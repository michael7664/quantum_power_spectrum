#pragma once
#include <vector>
#include <string>
#include "BackgroundCosmology.h"

struct Galaxy {
    double x, y, z;    // comoving Cartesian [Mpc/h]
    double weight;
};

class GalaxyCatalog {
public:
    std::vector<Galaxy> galaxies;
    std::string         source_name;
    std::string         coord_type;
    double              boxsize;
    double              redshift;
    bool                is_periodic;

    GalaxyCatalog() = default;

    // Load universal .bin written by tng_data_prep
    void load_bin(const std::string &filename,
                  const BackgroundCosmology &cosmo);

    // Load plain ASCII: x y z [weight]
    void load_ascii(const std::string &filename);

    size_t size() const { return galaxies.size(); }
    void   print_summary() const;

private:
    void convert_radecz(const BackgroundCosmology &cosmo);
};
