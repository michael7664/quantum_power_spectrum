#pragma once
#include "BackgroundCosmology.h"

class TransferFunction {
public:
    explicit TransferFunction(const BackgroundCosmology &cosmo);

    // Eisenstein & Hu (1998) full transfer function
    double T_EH(double k) const;       // k in h/Mpc

    // BBKS (1986) transfer function
    double T_BBKS(double k) const;

    // Primordial power spectrum P_prim(k) = A_s * (k/k_pivot)^(n_s-1)
    double P_primordial(double k) const;

    // Linear matter power spectrum P(k) = P_prim(k) * T(k)^2
    double P_linear(double k) const;

    double ns;          // spectral index
    double As;          // amplitude
    double k_pivot;     // pivot scale [h/Mpc]

private:
    const BackgroundCosmology &cosmo_;
    double z_eq_;       // matter-radiation equality redshift
    double k_eq_;       // equality scale [h/Mpc]
    double sound_horizon_;
    void   precompute();
};
