#include "LandySzalay.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>

LandySzalay::LandySzalay(double r_min, double r_max, int N_bins)
    : r_min_(r_min), r_max_(r_max), N_bins_(N_bins) {}

int LandySzalay::bin_index(double r) const {
    if (r < r_min_ || r >= r_max_) return -1;
    double log_r    = std::log(r);
    double log_rmin = std::log(r_min_);
    double log_rmax = std::log(r_max_);
    int    bin      = (int)((log_r - log_rmin)
                    / (log_rmax - log_rmin) * N_bins_);
    if (bin <  0)       return 0;
    if (bin >= N_bins_) return N_bins_ - 1;
    return bin;
}

GalaxyCatalog LandySzalay::make_randoms(
    const GalaxyCatalog &data,
    int N_randoms, int seed) const
{
    if (N_randoms < 0) N_randoms = 5 * (int)data.size();

    double xmin=1e30, xmax=-1e30;
    double ymin=1e30, ymax=-1e30;
    double zmin=1e30, zmax=-1e30;
    for (auto &g : data.galaxies) {
        xmin=std::min(xmin,g.x); xmax=std::max(xmax,g.x);
        ymin=std::min(ymin,g.y); ymax=std::max(ymax,g.y);
        zmin=std::min(zmin,g.z); zmax=std::max(zmax,g.z);
    }

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::uniform_real_distribution<double> uy(ymin, ymax);
    std::uniform_real_distribution<double> uz(zmin, zmax);

    GalaxyCatalog rcat;
    rcat.galaxies.reserve(N_randoms);
    rcat.source_name = "randoms";
    rcat.coord_type  = "cartesian_comoving";
    rcat.boxsize     = data.boxsize;
    rcat.redshift    = data.redshift;
    rcat.is_periodic = data.is_periodic;

    for (int i = 0; i < N_randoms; ++i)
        rcat.galaxies.push_back({ux(rng), uy(rng), uz(rng), 1.0});

    std::cout << "  Generated " << N_randoms << " randoms\n";
    return rcat;
}

std::vector<long long> LandySzalay::count_pairs(
    const GalaxyCatalog &A,
    const GalaxyCatalog &B,
    bool autocorr) const
{
    std::vector<long long> counts(N_bins_, 0);
    int NA = (int)A.size();
    int NB = (int)B.size();

#ifdef HAS_OPENMP
    #pragma omp parallel
    {
        std::vector<long long> local(N_bins_, 0);
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < NA; ++i) {
            int jstart = autocorr ? i+1 : 0;
            for (int j = jstart; j < NB; ++j) {
                double dx = A.galaxies[i].x - B.galaxies[j].x;
                double dy = A.galaxies[i].y - B.galaxies[j].y;
                double dz = A.galaxies[i].z - B.galaxies[j].z;
                double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
                int    b  = bin_index(r);
                if (b >= 0) local[b]++;
            }
        }
        #pragma omp critical
        for (int b = 0; b < N_bins_; ++b) counts[b] += local[b];
    }
#else
    for (int i = 0; i < NA; ++i) {
        int jstart = autocorr ? i+1 : 0;
        for (int j = jstart; j < NB; ++j) {
            double dx = A.galaxies[i].x - B.galaxies[j].x;
            double dy = A.galaxies[i].y - B.galaxies[j].y;
            double dz = A.galaxies[i].z - B.galaxies[j].z;
            double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
            int    b  = bin_index(r);
            if (b >= 0) counts[b]++;
        }
    }
#endif

    return counts;
}

std::vector<XiPoint> LandySzalay::compute(
    const GalaxyCatalog &data,
    const GalaxyCatalog &randoms) const
{
    long long ND = (long long)data.size();
    long long NR = (long long)randoms.size();

    // ── Pair count normalization factors ─────────────────────────────
    // DD: unique pairs from data    → ND*(ND-1)/2
    // DR: all cross pairs           → ND*NR
    // RR: unique pairs from randoms → NR*(NR-1)/2
    double n_DD = 0.5 * ND * (ND - 1);
    double n_DR = (double)ND * (double)NR;
    double n_RR = 0.5 * NR * (NR - 1);

    std::cout << "  ND=" << ND << "  NR=" << NR << "\n";
    std::cout << "  n_DD=" << n_DD
              << "  n_DR=" << n_DR
              << "  n_RR=" << n_RR << "\n";

    std::cout << "  Pair counting DD...\n";
    auto DD = count_pairs(data,    data,    true);
    std::cout << "  Pair counting DR...\n";
    auto DR = count_pairs(data,    randoms, false);
    std::cout << "  Pair counting RR...\n";
    auto RR = count_pairs(randoms, randoms, true);

    std::vector<XiPoint> result;
    result.reserve(N_bins_);

    double log_rmin = std::log(r_min_);
    double log_rmax = std::log(r_max_);

    for (int b = 0; b < N_bins_; ++b) {
        double log_r = log_rmin
                     + (b + 0.5) * (log_rmax - log_rmin) / N_bins_;
        double r     = std::exp(log_r);

        double dd = (double)DD[b];
        double dr = (double)DR[b];
        double rr = (double)RR[b];

        // ── Normalized pair counts ────────────────────────────────────
        double DD_n = (n_DD > 0) ? dd / n_DD : 0.0;
        double DR_n = (n_DR > 0) ? dr / n_DR : 0.0;
        double RR_n = (n_RR > 0) ? rr / n_RR : 0.0;

        // ── Landy-Szalay estimator ────────────────────────────────────
        // ξ(r) = (DD - 2·DR + RR) / RR   [normalized counts]
        double xi  = 0.0;
        double err = 0.0;
        if (RR_n > 0.0) {
            xi  = (DD_n - 2.0*DR_n + RR_n) / RR_n;
            // Poisson error estimate: σ = (1 + ξ) / sqrt(DD_raw)
            err = (dd > 0.0) ? (1.0 + xi) / std::sqrt(dd) : 0.0;
        }

        result.push_back({r, xi, err, DD[b], DR[b], RR[b]});
    }
    return result;
}

void LandySzalay::write(const std::vector<XiPoint> &xi,
                        const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# r[Mpc/h]  xi_LS(r)  error  DD  DR  RR\n";
    for (auto &p : xi)
        f << p.r   << "  " << p.xi    << "  " << p.error << "  "
          << p.DD  << "  " << p.DR    << "  " << p.RR    << "\n";
    std::cout << "  Written: " << filename << "\n";
}
