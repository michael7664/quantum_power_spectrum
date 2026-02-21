#pragma once
#include <vector>
#include <string>
#include "GalaxyCatalog.h"

struct XiPoint { double r, xi, error; long long DD, DR, RR; };

class LandySzalay {
public:
    LandySzalay(double r_min = 1.0,
                double r_max = 200.0,
                int    N_bins = 25);

    // Generate uniform random catalog in same volume
    GalaxyCatalog make_randoms(const GalaxyCatalog &data,
                               int    N_randoms = -1,
                               int    seed      = 42) const;

    // Compute ξ(r) using Landy-Szalay estimator
    std::vector<XiPoint> compute(const GalaxyCatalog &data,
                                 const GalaxyCatalog &randoms) const;

    static void write(const std::vector<XiPoint> &xi,
                      const std::string &filename);

private:
    double r_min_, r_max_;
    int    N_bins_;

    int bin_index(double r) const;
    std::vector<long long> count_pairs(
        const GalaxyCatalog &A,
        const GalaxyCatalog &B,
        bool autocorr) const;
};