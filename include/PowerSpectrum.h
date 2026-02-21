#pragma once
#include <vector>
#include <string>
#include "BackgroundCosmology.h"
#include "TransferFunction.h"

struct PKPoint { double k, Pk; };

class PowerSpectrum {
public:
    PowerSpectrum(const BackgroundCosmology &cosmo,
                  const TransferFunction    &tf);

    // Compute theory P(k) on a log-spaced grid
    // b_halo: linear halo bias (set to 1.0 for matter P(k))
    std::vector<PKPoint> compute(
        double k_min  = 1e-3,
        double k_max  = 10.0,
        int    N_bins = 200,
        double b_halo = 1.0) const;

    // Estimate linear bias from catalog mass
    // Uses Tinker et al. 2010 fitting formula
    static double bias_from_mass(double M200c_msun, double z = 0.0);

    // Correlation function ξ(r) via direct Fourier transform
    // handles zero-crossings correctly (no log scale artefacts)
    std::vector<std::pair<double,double>> xi_theory(
        double r_min  = 1.0,
        double r_max  = 250.0,
        int    N_r    = 100,
        double b_halo = 1.0) const;

    // Write P(k) to file
    static void write(const std::vector<PKPoint> &pk,
                      const std::string &filename);

    // Write ξ(r) to file
    static void write_xi(
        const std::vector<std::pair<double,double>> &xi,
        const std::string &filename);

private:
    const BackgroundCosmology &cosmo_;
    const TransferFunction    &tf_;

    double sigma8_norm() const;   // returns normalisation factor
};