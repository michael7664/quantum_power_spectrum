#pragma once
#include <vector>
#include <complex>
#include <string>
#include "BackgroundCosmology.h"

struct QFTCorrection {
    double k;
    double delta_P;       // quantum correction to P(k)
    double P_classical;
    double P_quantum;
};

class QuantumQFT {
public:
    explicit QuantumQFT(const BackgroundCosmology &cosmo);

    // One-loop quantum correction to P(k) from scalar field vacuum fluctuations
    double quantum_correction(double k, double Pk_classical) const;

    // Apply QFT corrections to a full P(k) vector
    std::vector<QFTCorrection> apply(
        const std::vector<std::pair<double,double>> &pk_classical) const;

    // Quantum-corrected ξ(r) from corrected P(k)
    std::vector<std::pair<double,double>> xi_quantum(
        const std::vector<QFTCorrection> &qft,
        double r_min = 1.0,
        double r_max = 200.0,
        int    N_r   = 60) const;

    static void write(const std::vector<QFTCorrection> &qft,
                      const std::string &filename);

    static void write_xi(
        const std::vector<std::pair<double,double>> &xi,
        const std::string &filename);

    // QFT parameters
    double lambda_;       // quartic coupling
    double m_phi_;        // scalar field mass [h/Mpc]
    double hbar_;         // effective ℏ (set to 1 in natural units)

private:
    const BackgroundCosmology &cosmo_;
    double loop_integral(double k) const;
};