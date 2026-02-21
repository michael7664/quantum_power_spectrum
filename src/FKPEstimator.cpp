#include "FKPEstimator.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <complex>

#ifdef HAS_FFTW3
  #include <fftw3.h>
#endif

static constexpr double PI = M_PI;

FKPEstimator::FKPEstimator(double P0, double k_min, double k_max,
                           int N_bins, int N_grid)
    : P0_(P0), k_min_(k_min), k_max_(k_max),
      N_bins_(N_bins), N_grid_(N_grid) {}

double FKPEstimator::fkp_weight(double nbar) const {
    return 1.0 / (1.0 + nbar * P0_);
}

// ── CIC (Cloud-In-Cell) mass assignment ───────────────────────────────
void FKPEstimator::assign_cic(
    const GalaxyCatalog &cat,
    std::vector<double>  &grid,
    int Ng,
    double x0, double y0, double z0,
    double Lx, double Ly, double Lz) const
{
    double dx = Lx / Ng;
    double dy = Ly / Ng;
    double dz = Lz / Ng;

    for (auto &g : cat.galaxies) {
        double fx = (g.x - x0) / dx;
        double fy = (g.y - y0) / dy;
        double fz = (g.z - z0) / dz;

        int ix = (int)std::floor(fx);
        int iy = (int)std::floor(fy);
        int iz = (int)std::floor(fz);

        double wx1 = fx - ix;  double wx0 = 1.0 - wx1;
        double wy1 = fy - iy;  double wy0 = 1.0 - wy1;
        double wz1 = fz - iz;  double wz0 = 1.0 - wz1;

        for (int ddx = 0; ddx <= 1; ++ddx)
        for (int ddy = 0; ddy <= 1; ++ddy)
        for (int ddz = 0; ddz <= 1; ++ddz) {
            int cx = std::max(0, std::min(ix + ddx, Ng - 1));
            int cy = std::max(0, std::min(iy + ddy, Ng - 1));
            int cz = std::max(0, std::min(iz + ddz, Ng - 1));

            double w = (ddx==0 ? wx0 : wx1)
                     * (ddy==0 ? wy0 : wy1)
                     * (ddz==0 ? wz0 : wz1)
                     * g.weight;

            grid[cx * Ng*Ng + cy * Ng + cz] += w;
        }
    }
}

// ── Spherical binning of |F(k)|² ─────────────────────────────────────
std::vector<PKData> FKPEstimator::bin_power(
    const std::vector<std::complex<double>> &Fk,
    int Ng, double Lbox, double shot_noise) const
{
    double dk_fund  = 2.0 * PI / Lbox;
    double Vcell    = std::pow(Lbox / Ng, 3.0);
    double norm     = Vcell * Vcell / (Lbox * Lbox * Lbox);

    std::vector<double> Pk_sum (N_bins_, 0.0);
    std::vector<int>    Nmodes (N_bins_, 0);

    double log_kmin = std::log(k_min_);
    double log_kmax = std::log(k_max_);
    double dlogk    = (log_kmax - log_kmin) / N_bins_;

    int Ng2 = Ng / 2;

    for (int ix = 0; ix < Ng; ++ix)
    for (int iy = 0; iy < Ng; ++iy)
    for (int iz = 0; iz < Ng; ++iz) {
        int fx = (ix <= Ng2) ? ix : ix - Ng;
        int fy = (iy <= Ng2) ? iy : iy - Ng;
        int fz = (iz <= Ng2) ? iz : iz - Ng;

        double k = dk_fund * std::sqrt(
            (double)(fx*fx + fy*fy + fz*fz));

        if (k < k_min_ || k >= k_max_) continue;

        int bin = (int)((std::log(k) - log_kmin) / dlogk);
        if (bin < 0 || bin >= N_bins_) continue;

        int    idx   = ix * Ng*Ng + iy * Ng + iz;
        double power = std::norm(Fk[idx]) * norm;

        Pk_sum[bin] += power;
        Nmodes[bin] += 1;
    }

    std::vector<PKData> result;
    result.reserve(N_bins_);

    for (int b = 0; b < N_bins_; ++b) {
        double k   = std::exp(log_kmin + (b + 0.5) * dlogk);
        double Pk  = (Nmodes[b] > 0)
                   ? Pk_sum[b] / Nmodes[b] - shot_noise
                   : 0.0;
        double Pk_sub = std::max(0.0, Pk);
        result.push_back({k, Pk, Pk_sub, Nmodes[b]});
    }
    return result;
}

// ── Full FFT estimator ────────────────────────────────────────────────
std::vector<PKData> FKPEstimator::compute(
    const GalaxyCatalog       &data,
    const GalaxyCatalog       &randoms,
    const BackgroundCosmology &cosmo) const
{
#ifndef HAS_FFTW3
    std::cout << "  [FKP] WARNING: FFTW3 not compiled in"
              << " — using direct estimator\n";
    return compute_direct(data, randoms, cosmo);
#else
    std::cout << "  [FKP] Using FFTW3 grid estimator (Ng="
              << N_grid_ << ")\n";

    int    ND    = (int)data.size();
    int    NR    = (int)randoms.size();
    double alpha = (double)ND / (double)NR;
    int    Ng    = N_grid_;
    int    Ng3   = Ng * Ng * Ng;

    // ── Bounding box with 5% padding ─────────────────────────────
    double xmin=1e30, xmax=-1e30;
    double ymin=1e30, ymax=-1e30;
    double zmin=1e30, zmax=-1e30;
    for (auto &g : data.galaxies) {
        xmin=std::min(xmin,g.x); xmax=std::max(xmax,g.x);
        ymin=std::min(ymin,g.y); ymax=std::max(ymax,g.y);
        zmin=std::min(zmin,g.z); zmax=std::max(zmax,g.z);
    }
    double Lx   = (xmax - xmin) * 1.05;
    double Ly   = (ymax - ymin) * 1.05;
    double Lz   = (zmax - zmin) * 1.05;
    double Lbox = std::max({Lx, Ly, Lz});
    double x0   = 0.5*(xmin+xmax) - 0.5*Lbox;
    double y0   = 0.5*(ymin+ymax) - 0.5*Lbox;
    double z0   = 0.5*(zmin+zmax) - 0.5*Lbox;

    double V_survey = Lbox * Lbox * Lbox;
    double nbar     = ND / V_survey;
    double w_fkp    = fkp_weight(nbar);

    // Shot noise: P_shot = (1 + alpha) / nbar
    double shot = (1.0 + alpha) / nbar;

    std::cout << "  [FKP] Lbox     = " << Lbox      << " Mpc/h\n";
    std::cout << "  [FKP] nbar     = " << nbar      << " (Mpc/h)^-3\n";
    std::cout << "  [FKP] w_FKP    = " << w_fkp     << "\n";
    std::cout << "  [FKP] shot     = " << shot       << " (Mpc/h)^3\n";
    std::cout << "  [FKP] alpha    = " << alpha      << "\n";

    // ── Grid galaxies and randoms ─────────────────────────────────
    std::vector<double> grid_g(Ng3, 0.0);
    std::vector<double> grid_r(Ng3, 0.0);
    assign_cic(data,    grid_g, Ng, x0, y0, z0, Lbox, Lbox, Lbox);
    assign_cic(randoms, grid_r, Ng, x0, y0, z0, Lbox, Lbox, Lbox);

    // ── F(r) = w_FKP * [n_g(r) - alpha * n_r(r)] ─────────────────
    std::vector<double> Fr(Ng3);
    for (int i = 0; i < Ng3; ++i)
        Fr[i] = w_fkp * (grid_g[i] - alpha * grid_r[i]);

    // ── FFT: F(r) → F(k) ─────────────────────────────────────────
    fftw_complex *in  = fftw_alloc_complex(Ng3);
    fftw_complex *out = fftw_alloc_complex(Ng3);

    for (int i = 0; i < Ng3; ++i) {
        in[i][0] = Fr[i];
        in[i][1] = 0.0;
    }

    fftw_plan plan = fftw_plan_dft_3d(
        Ng, Ng, Ng, in, out,
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Copy output
    std::vector<std::complex<double>> Fk(Ng3);
    for (int i = 0; i < Ng3; ++i)
        Fk[i] = {out[i][0], out[i][1]};

    fftw_free(in);
    fftw_free(out);

    // ── Spherically average into k-bins ──────────────────────────
    return bin_power(Fk, Ng, Lbox, shot);
#endif
}

// ── Direct estimator (no FFT fallback) ───────────────────────────────
std::vector<PKData> FKPEstimator::compute_direct(
    const GalaxyCatalog       &data,
    const GalaxyCatalog       &randoms,
    const BackgroundCosmology &cosmo) const
{
    int    ND    = (int)data.size();
    int    NR    = (int)randoms.size();
    double alpha = (double)ND / (double)NR;

    double xmin=1e30, xmax=-1e30;
    double ymin=1e30, ymax=-1e30;
    double zmin=1e30, zmax=-1e30;
    for (auto &g : data.galaxies) {
        xmin=std::min(xmin,g.x); xmax=std::max(xmax,g.x);
        ymin=std::min(ymin,g.y); ymax=std::max(ymax,g.y);
        zmin=std::min(zmin,g.z); zmax=std::max(zmax,g.z);
    }
    double V_survey = (xmax-xmin) * (ymax-ymin) * (zmax-zmin);
    double nbar     = ND / V_survey;
    double shot     = (1.0 + alpha) / nbar;

    std::cout << "  [FKP direct] V=" << V_survey
              << "  nbar=" << nbar
              << "  shot=" << shot << "\n";

    std::vector<PKData> result;
    result.reserve(N_bins_);

    double log_kmin = std::log(k_min_);
    double log_kmax = std::log(k_max_);
    double dlogk    = (log_kmax - log_kmin) / N_bins_;

    // Cap at 500 galaxies for speed
    int Nmax = std::min(ND, 500);

    for (int b = 0; b < N_bins_; ++b) {
        double k    = std::exp(log_kmin + (b + 0.5) * dlogk);
        double k_lo = std::exp(log_kmin +  b        * dlogk);
        double k_hi = std::exp(log_kmin + (b + 1)   * dlogk);

        double V_shell = 4.0/3.0 * PI
                       * (k_hi*k_hi*k_hi - k_lo*k_lo*k_lo);
        int N_modes = std::max(1,
            (int)(V_shell * V_survey / (8.0 * PI*PI*PI)));

        double Pk_sum = 0.0;
        int    N_pairs = 0;

        for (int i = 0; i < Nmax; ++i)
        for (int j = i+1; j < Nmax; ++j) {
            double dx = data.galaxies[i].x - data.galaxies[j].x;
            double dy = data.galaxies[i].y - data.galaxies[j].y;
            double dz = data.galaxies[i].z - data.galaxies[j].z;
            double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-6) continue;
            double kr   = k * r;
            double sinc = (kr > 1e-6) ? std::sin(kr)/kr
                                      : 1.0 - kr*kr/6.0;
            Pk_sum += sinc;
            N_pairs++;
        }

        double Pk_raw = (N_pairs > 0)
                      ? 2.0 * V_survey * Pk_sum
                        / ((double)N_pairs * (double)ND)
                      : 0.0;
        double Pk_sub = std::max(0.0, Pk_raw - shot);
        result.push_back({k, Pk_raw, Pk_sub, N_modes});
    }
    return result;
}

// ── File I/O ──────────────────────────────────────────────────────────
void FKPEstimator::write(const std::vector<PKData> &pk,
                         const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# k[h/Mpc]  P(k)  P(k)_shot_subtracted  N_modes\n";
    for (auto &p : pk)
        f << p.k << "  " << p.Pk << "  "
          << p.Pk_shot_subtracted << "  " << p.N_modes << "\n";
    std::cout << "  Written: " << filename << "\n";
}
