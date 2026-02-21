#pragma once
#include <vector>
#include <string>
#include <complex>
#include "GalaxyCatalog.h"
#include "BackgroundCosmology.h"

struct PKData {
    double k;
    double Pk;                   // raw estimated P(k)
    double Pk_shot_subtracted;   // shot-noise subtracted
    int    N_modes;              // number of Fourier modes in bin
};

class FKPEstimator {
public:
    FKPEstimator(double P0     = 1e4,    // FKP reference power [(Mpc/h)³]
                 double k_min  = 0.005,
                 double k_max  = 0.5,
                 int    N_bins = 30,
                 int    N_grid = 64);    // FFT grid size per side

    // Full FFT-based P(k) estimator (CIC mass assignment + FFT)
    std::vector<PKData> compute(
        const GalaxyCatalog       &data,
        const GalaxyCatalog       &randoms,
        const BackgroundCosmology &cosmo) const;

    // Simple direct estimator (fallback, no FFT required)
    std::vector<PKData> compute_direct(
        const GalaxyCatalog       &data,
        const GalaxyCatalog       &randoms,
        const BackgroundCosmology &cosmo) const;

    static void write(const std::vector<PKData> &pk,
                      const std::string &filename);

private:
    double P0_;
    double k_min_;
    double k_max_;
    int    N_bins_;
    int    N_grid_;

    double fkp_weight(double nbar) const;

    // CIC (Cloud-In-Cell) mass assignment to 3D grid
    void assign_cic(
        const GalaxyCatalog &cat,
        std::vector<double>  &grid,
        int Ng,
        double x0, double y0, double z0,
        double Lx, double Ly, double Lz) const;

    // Spherically average |F(k)|² into k-bins
    std::vector<PKData> bin_power(
        const std::vector<std::complex<double>> &Fk,
        int Ng, double Lbox,
        double shot_noise) const;
};
