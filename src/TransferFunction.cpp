#include "TransferFunction.h"
#include <cmath>

static constexpr double PI = M_PI;

TransferFunction::TransferFunction(const BackgroundCosmology &cosmo)
    : cosmo_(cosmo), ns(0.9649), As(2.1e-9), k_pivot(0.05)
{
    precompute();
}

void TransferFunction::precompute() {
    double h  = cosmo_.h;
    double Oc = cosmo_.OmegaM;
    double Ob = 0.0486;               // baryon fraction (Planck 2018)
    double T  = cosmo_.Tcmb / 2.7;

    // Matter-radiation equality
    z_eq_ = 2.50e4 * Oc * h * h * std::pow(T, -4);
    k_eq_ = 7.46e-2 * Oc * h * h / (T * T);     // h/Mpc

    // Sound horizon (Eisenstein & Hu eq. 6)
    double z_d   = 1291.0 * std::pow(Oc*h*h, 0.251)
                 / (1.0 + 0.659 * std::pow(Oc*h*h, 0.828))
                 * (1.0 + std::pow(Ob*h*h, 0.2) * 0.828);
    double R_eq  = 31.5e3 * Ob*h*h * std::pow(T,-4) * (1000.0/z_eq_);
    double R_d   = 31.5e3 * Ob*h*h * std::pow(T,-4) * (1000.0/z_d);
    sound_horizon_ = 2.0 / (3.0 * k_eq_)
                   * std::sqrt(6.0 / R_eq)
                   * std::log((std::sqrt(1+R_d) + std::sqrt(R_d+R_eq))
                              / (1.0 + std::sqrt(R_eq)));
}

double TransferFunction::T_EH(double k) const {
    // Eisenstein & Hu (1998) zero-baryon approximation
    double h   = cosmo_.h;
    double Om  = cosmo_.OmegaM;
    double T27 = cosmo_.Tcmb / 2.7;
    double q   = k / (13.41 * k_eq_);
    double C0  = 14.2 + 731.0 / (1.0 + 62.5 * q);
    double T0  = std::log(M_E + 1.8*q) / (std::log(M_E + 1.8*q) + C0*q*q);

    // Silk damping
    double k_silk = 1.6 * std::pow(Om*h*h, 0.52)
                        * std::pow(0.0486*h*h, 0.09)  // Ob h^2
                        * std::pow(T27, -0.95);        // h/Mpc
    double f      = 1.0 / (1.0 + std::pow(k * sound_horizon_ / 5.4, 4));
    double C1     = 14.2 / 1.0 + 386.0 / (1.0 + 69.9 * std::pow(q, 1.08));
    double T1     = std::log(M_E + 1.8*q) / (std::log(M_E + 1.8*q) + C1*q*q);
    return f * T1 + (1.0 - f) * T0;
}

double TransferFunction::T_BBKS(double k) const {
    double h  = cosmo_.h;
    double Om = cosmo_.OmegaM;
    double q  = k / (Om * h * h) * std::exp(cosmo_.OmegaR + 0.0486);
    return std::log(1.0 + 2.34*q) / (2.34*q)
         * std::pow(1.0 + 3.89*q
           + std::pow(16.1*q,2)
           + std::pow(5.46*q,3)
           + std::pow(6.71*q,4), -0.25);
}

double TransferFunction::P_primordial(double k) const {
    return As * std::pow(k / k_pivot, ns - 1.0);
}

double TransferFunction::P_linear(double k) const {
    double T = T_EH(k);
    return P_primordial(k) * T * T * k;   // k factor: Harrison-Zel'dovich
}