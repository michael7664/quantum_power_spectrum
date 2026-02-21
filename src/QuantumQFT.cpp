#include "QuantumQFT.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>

static constexpr double PI = M_PI;

QuantumQFT::QuantumQFT(const BackgroundCosmology &cosmo)
    : cosmo_(cosmo),
      lambda_(0.1),
      m_phi_(0.01),
      hbar_(1.0)
{}

// ── One-loop scalar field correction ─────────────────────────────────
// Pauli-Villars regulated result:
//   ΔP(k)/P(k) = λ/(16π²) · ln(1 + Λ²/(m² + k²))
// Dimensionless by construction, ~0.67% at k→0, decays to 0 at k→Λ
double QuantumQFT::loop_integral(double k) const {
    double Lambda = 2.0;    // UV cutoff [h/Mpc]
    double m2     = m_phi_ * m_phi_;
    double k2     = k * k;
    double Lam2   = Lambda * Lambda;
    return std::log(1.0 + Lam2 / (m2 + k2));
}

double QuantumQFT::quantum_correction(double k,
                                      double Pk_classical) const {
    double prefactor = lambda_ * hbar_ / (16.0 * PI * PI);
    return prefactor * loop_integral(k) * Pk_classical;
}

// ── Apply to full P(k) vector ─────────────────────────────────────────
std::vector<QFTCorrection> QuantumQFT::apply(
    const std::vector<std::pair<double,double>> &pk_classical) const
{
    std::vector<QFTCorrection> result;
    result.reserve(pk_classical.size());
    for (auto &[k, Pk] : pk_classical) {
        double dP = quantum_correction(k, Pk);
        result.push_back({k, dP, Pk, Pk + dP});
    }
    return result;
}

// ── Quantum ξ(r) via Fourier transform of P_quantum(k) ───────────────
std::vector<std::pair<double,double>> QuantumQFT::xi_quantum(
    const std::vector<QFTCorrection> &qft,
    double r_min, double r_max, int N_r) const
{
    std::vector<std::pair<double,double>> result;
    result.reserve(N_r);

    double dr = (r_max - r_min) / (N_r - 1);
    int    N  = (int)qft.size();

    for (int ir = 0; ir < N_r; ++ir) {
        double r  = r_min + ir * dr;
        double xi = 0.0;

        for (int ik = 0; ik < N - 1; ++ik) {
            double k1 = qft[ik  ].k;
            double k2 = qft[ik+1].k;
            double P1 = qft[ik  ].P_quantum;
            double P2 = qft[ik+1].P_quantum;

            if (P1 <= 0.0 || P2 <= 0.0) continue;

            double dk = k2 - k1;
            double km = 0.5 * (k1 + k2);
            double Pm = 0.5 * (P1 + P2);
            double kr = km * r;

            // sinc: well-behaved at kr → 0
            double sinc = (kr > 1e-6)
                        ? std::sin(kr) / kr
                        : 1.0 - kr*kr/6.0 + kr*kr*kr*kr/120.0;

            // Hann taper: suppresses Gibbs ringing from sharp k cutoff
            constexpr double k_taper_lo = 5.0;
            constexpr double k_taper_hi = 10.0;
            double taper = (km < k_taper_lo) ? 1.0
                         : (km > k_taper_hi) ? 0.0
                         : 0.5 * (1.0 + std::cos(
                               M_PI * (km - k_taper_lo)
                                     / (k_taper_hi - k_taper_lo)));

            xi += Pm * sinc * km * km * dk * taper;
        }
        xi /= 2.0 * PI * PI;
        result.push_back({r, xi});
    }
    return result;
}

// ── File I/O ──────────────────────────────────────────────────────────
void QuantumQFT::write(const std::vector<QFTCorrection> &qft,
                       const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# k[h/Mpc]  delta_P  P_classical  P_quantum\n";
    for (auto &q : qft)
        f << q.k          << "  " << q.delta_P     << "  "
          << q.P_classical << "  " << q.P_quantum   << "\n";
    std::cout << "  Written: " << filename << "\n";
}

void QuantumQFT::write_xi(
    const std::vector<std::pair<double,double>> &xi,
    const std::string &filename)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot write: " + filename);
    f << "# r[Mpc/h]   xi_quantum(r)\n";
    for (auto &p : xi)
        f << p.first << "   " << p.second << "\n";
    std::cout << "  Written: " << filename << "\n";
}
