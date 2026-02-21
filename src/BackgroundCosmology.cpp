#include "BackgroundCosmology.h"
#include <iostream>
#include <stdexcept>
#include <cmath>

static constexpr double PI     = M_PI;
static constexpr int    N_INTE = 1000;   // integration steps

BackgroundCosmology::BackgroundCosmology(
    double H0_, double OmegaM_, double OmegaL_,
    double OmegaR_, double Tcmb_)
    : H0(H0_), h(H0_ / 100.0),
      OmegaM(OmegaM_), OmegaL(OmegaL_),
      OmegaR(OmegaR_), Tcmb(Tcmb_)
{
    OmegaK = 1.0 - OmegaM - OmegaL - OmegaR;
}

double BackgroundCosmology::H(double z) const {
    double E2 = OmegaM * std::pow(1+z, 3)
              + OmegaR * std::pow(1+z, 4)
              + OmegaL
              + OmegaK * std::pow(1+z, 2);
    return H0 * std::sqrt(E2);
}

double BackgroundCosmology::integrate_chi(double z) const {
    // Trapezoidal rule:  chi = c/H0 * integral_0^z dz'/E(z')
    // Returns chi in Mpc/h  (c/H0 in units where c=1 → c=2997.92 Mpc/h)
    constexpr double c_over_H0 = 2997.92;   // c [km/s] / 100 [km/s/Mpc]
    if (z <= 0.0) return 0.0;

    double dz   = z / N_INTE;
    double sum  = 0.0;
    double zprev = 0.0, fprev = 1.0 / (H(0) / H0);

    for (int i = 1; i <= N_INTE; ++i) {
        double zi = i * dz;
        double fi = H0 / H(zi);
        sum   += 0.5 * (fprev + fi) * dz;
        zprev  = zi;
        fprev  = fi;
    }
    return c_over_H0 * sum;
}

double BackgroundCosmology::chi(double z) const {
    return integrate_chi(z);
}

double BackgroundCosmology::dA(double z) const {
    return chi(z) / (1.0 + z);
}

double BackgroundCosmology::dL(double z) const {
    return chi(z) * (1.0 + z);
}

double BackgroundCosmology::integrate_growth(double z) const {
    // D(z) ∝ H(z) * integral_z^inf dz'/(H(z')^3 * (1+z')^... )
    // Simplified Carroll et al. 1992 approximation
    double Om = OmegaM * std::pow(1+z, 3);
    double OL = OmegaL;
    double E2 = Om + OL + OmegaK * std::pow(1+z, 2)
                        + OmegaR * std::pow(1+z, 4);
    double E  = std::sqrt(E2 / (H0 * H0));
    // Integral approximation (Heath 1977)
    constexpr double c_H0 = 2997.92;
    double sum  = 0.0;
    int    Nint = 2000;
    double zmax = 1000.0;
    double dz   = (zmax - z) / Nint;
    double zp   = z;
    double fp   = 1.0 / std::pow(H(zp)/H0 * (1+zp), 3.0);

    for (int i = 1; i <= Nint; ++i) {
        zp      = z + i * dz;
        double fn = 1.0 / std::pow(H(zp)/H0 * (1+zp), 3.0);
        sum    += 0.5 * (fp + fn) * dz;
        fp      = fn;
    }
    return (H(z) / H0) * sum;
}

double BackgroundCosmology::growthFactor(double z) const {
    double D0 = integrate_growth(0.0);
    double Dz = integrate_growth(z);
    return (D0 > 0) ? Dz / D0 : 1.0;
}

void BackgroundCosmology::radecz_to_xyz(
    double ra_deg, double dec_deg, double z,
    double &x, double &y, double &z_out) const
{
    constexpr double DEG2RAD = PI / 180.0;
    double r   = chi(z);
    double ra  = ra_deg  * DEG2RAD;
    double dec = dec_deg * DEG2RAD;
    x     = r * std::cos(dec) * std::cos(ra);
    y     = r * std::cos(dec) * std::sin(ra);
    z_out = r * std::sin(dec);
}

void BackgroundCosmology::print() const {
    std::cout << "\n── BackgroundCosmology ──────────────────\n";
    std::cout << "  H0     = " << H0     << " km/s/Mpc\n";
    std::cout << "  OmegaM = " << OmegaM << "\n";
    std::cout << "  OmegaL = " << OmegaL << "\n";
    std::cout << "  OmegaR = " << OmegaR << "\n";
    std::cout << "  OmegaK = " << OmegaK << "\n";
    std::cout << "  chi(z=0.5) = " << chi(0.5) << " Mpc/h\n";
    std::cout << "  chi(z=1.0) = " << chi(1.0) << " Mpc/h\n";
    std::cout << "  D(z=1)/D(0)= " << growthFactor(1.0) << "\n";
    std::cout << "─────────────────────────────────────────\n";
}
