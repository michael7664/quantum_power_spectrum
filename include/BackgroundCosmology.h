#pragma once
#include <cmath>
#include <vector>
#include <string>

class BackgroundCosmology {
public:
    // Cosmological parameters (Planck 2018)
    double H0;        // Hubble constant [km/s/Mpc]
    double h;         // H0 / 100
    double OmegaM;    // matter density
    double OmegaL;    // dark energy density
    double OmegaR;    // radiation density
    double OmegaK;    // curvature (derived)
    double Tcmb;      // CMB temperature [K]

    BackgroundCosmology(
        double H0    = 67.4,
        double OmegaM = 0.315,
        double OmegaL = 0.685,
        double OmegaR = 9.0e-5,
        double Tcmb   = 2.7255);

    // Hubble parameter at redshift z [km/s/Mpc]
    double H(double z) const;

    // Comoving distance to redshift z [Mpc/h]
    double chi(double z) const;

    // Angular diameter distance [Mpc/h]
    double dA(double z) const;

    // Luminosity distance [Mpc/h]
    double dL(double z) const;

    // Growth factor D(z) normalised to D(0)=1
    double growthFactor(double z) const;

    // Convert (RA, DEC, z) → comoving Cartesian (x,y,z) [Mpc/h]
    void radecz_to_xyz(double ra_deg, double dec_deg, double z,
                       double &x, double &y, double &z_out) const;

    void print() const;

private:
    // Numerical integration helpers
    double integrate_chi(double z) const;
    double integrate_growth(double z) const;
};