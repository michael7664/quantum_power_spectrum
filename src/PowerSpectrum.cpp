#include "PowerSpectrum.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>

static constexpr double PI = M_PI;

PowerSpectrum::PowerSpectrum(const BackgroundCosmology &cosmo,
                             const TransferFunction    &tf)
    : cosmo_(cosmo), tf_(tf) {}

// ── sigma8 normalisation ──────────────────────────────────────────────
double PowerSpectrum::sigma8_norm() const {
    double sigma8_target = 0.811;
    double R8     = 8.0;
    int    Nint   = 2000;
    double lkmin  = std::log(1e-4);
    double lkmax  = std::log(50.0);
    double dlk    = (lkmax - lkmin) / Nint;
    double sum    = 0.0;

    for (int i = 0; i < Nint; ++i) {
        double k  = std::exp(lkmin + (i + 0.5) * dlk);
        double kR = k * R8;
        double W  = 3.0 * (std::sin(kR) - kR * std::cos(kR))
                  / (kR * kR * kR);
        double Pk = tf_.P_linear(k);
        sum      += Pk * W * W * k * k * k * dlk / (2.0 * PI * PI);
    }
    double sigma8_raw = std::sqrt(sum);
    return (sigma8_raw > 0)
         ? (sigma8_target * sigma8_target) / (sigma8_raw * sigma8_raw)
         : 1.0;
}

// ── Tinker et al. 2010 linear bias ───────────────────────────────────
double PowerSpectrum::bias_from_mass(double M200c_msun, double z) {
    constexpr double delta_c = 1.686;
    constexpr double q       = 0.707;
    constexpr double p       = 0.3;

    double M_star = 3.0e12;
    double alpha  = 0.3;
    double sigma  = std::pow(M200c_msun / M_star, -alpha);
    double nu     = delta_c / sigma;

    double b = 1.0
             + (q * nu * nu - 1.0) / delta_c
             + 2.0 * p / delta_c
               / (1.0 + std::pow(q * nu * nu, p));
    return std::max(1.0, b);
}

// ── P(k) with bias ───────────────────────────────────────────────────
std::vector<PKPoint> PowerSpectrum::compute(
    double k_min, double k_max, int N_bins, double b_halo) const
{
    std::vector<PKPoint> result;
    result.reserve(N_bins);

    double norm     = sigma8_norm();
    double b2       = b_halo * b_halo;
    double log_kmin = std::log(k_min);
    double log_kmax = std::log(k_max);
    double dlogk    = (log_kmax - log_kmin) / (N_bins - 1);

    for (int i = 0; i < N_bins; ++i) {
        double k  = std::exp(log_kmin + i * dlogk);
        double Pk = tf_.P_linear(k) * norm * b2;
        result.push_back({k, Pk});
    }
    return result;
}

// ── ξ(r) via direct Fourier transform ────────────────────────────────
std::vector<std::pair<double,double>> PowerSpectrum::xi_theory(
    double r_min, double r_max, int N_r, double b_halo) const
{
    std::vector<std::pair<double,double>> result;
    result.reserve(N_r);

    double norm  = sigma8_norm();
    double b2    = b_halo * b_halo;

    // High-resolution unbiased P(k) grid (bias applied via b2 below)
    // Extended to k=50 h/Mpc with 8000 points for smooth xi(r)
    auto pk_vec = compute(1e-4, 10.0, 4000, 1.0);

    double dr = (r_max - r_min) / (N_r - 1);
    int    N  = (int)pk_vec.size();

    for (int ir = 0; ir < N_r; ++ir) {
        double r  = r_min + ir * dr;
        double xi = 0.0;

        for (int ik = 0; ik < N - 1; ++ik) {
            double k1 = pk_vec[ik  ].k;
            double k2 = pk_vec[ik+1].k;
            double P1 = pk_vec[ik  ].Pk;
            double P2 = pk_vec[ik+1].Pk;
            double dk = k2 - k1;
            double km = 0.5 * (k1 + k2);
            double Pm = 0.5 * (P1 + P2);
            double kr = km * r;

            // sinc: well-behaved at kr → 0
            double sinc = (kr > 1e-6)
                        ? std::sin(kr) / kr
                        : 1.0 - kr*kr/6.0 + kr*kr*kr*kr/120.0;

            // Hann taper: suppresses Gibbs ringing from sharp k cutoff
            // Smoothly damps modes between 5 and 50 h/Mpc to zero
            constexpr double k_taper_lo = 1.0;
            constexpr double k_taper_hi = 10.0;
            double taper = (km < k_taper_lo) ? 1.0
                         : (km > k_taper_hi) ? 0.0
                         : 0.5 * (1.0 + std::cos(
                               M_PI * (km - k_taper_lo)
                                     / (k_taper_hi - k_taper_lo)));

            xi += Pm * sinc * km * km * dk * taper;
        }
        xi *= b2 / (2.0 * PI * PI);
        result.push_back({r, xi});
    }
    return result;
}

// ── File I/O ──────────────────────────────────────────────────────────
void PowerSpectrum::write(const std::vector<PKPoint> &pk,
                          const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# k[h/Mpc]   P(k)[(Mpc/h)^3]\n";
    for (auto &p : pk)
        f << p.k << "   " << p.Pk << "\n";
    std::cout << "  Written: " << filename
              << "  (" << pk.size() << " points)\n";
}

void PowerSpectrum::write_xi(
    const std::vector<std::pair<double,double>> &xi,
    const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# r[Mpc/h]   xi(r)\n";
    for (auto &p : xi)
        f << p.first << "   " << p.second << "\n";
    std::cout << "  Written: " << filename
              << "  (" << xi.size() << " points)\n";
}