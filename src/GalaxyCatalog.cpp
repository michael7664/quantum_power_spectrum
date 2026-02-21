#include "GalaxyCatalog.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <cmath>

static constexpr int HEADER_SIZE  = 256;
static constexpr int RECORD_SIZE  = 32;   // 4 × double
static const     char MAGIC[8]    = {'Q','P','S','C','_','v','1','\0'};

void GalaxyCatalog::load_bin(const std::string &filename,
                              const BackgroundCosmology &cosmo)
{
    std::ifstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open: " + filename);

    // ── Read header ──────────────────────────────────────────────────
    char header[HEADER_SIZE] = {};
    f.read(header, HEADER_SIZE);

    // Verify magic
    if (std::memcmp(header, MAGIC, 8) != 0)
        throw std::runtime_error("Bad magic in: " + filename);

    int64_t N_gal;
    double  bsize, z_val;
    uint8_t is_per, has_w;

    std::memcpy(&N_gal, header +  8, 8);
    std::memcpy(&bsize, header + 16, 8);
    std::memcpy(&z_val, header + 24, 8);
    std::memcpy(&is_per,header + 32, 1);
    std::memcpy(&has_w, header + 33, 1);

    char sname[65] = {};   std::memcpy(sname, header + 34, 64);
    char ctype[65] = {};   std::memcpy(ctype, header + 98, 64);

    boxsize      = bsize;
    redshift     = z_val;
    is_periodic  = (bool)is_per;
    source_name  = std::string(sname);
    coord_type   = std::string(ctype);

    std::cout << "\n── Loading catalog ──────────────────────────\n";
    std::cout << "  File       : " << filename   << "\n";
    std::cout << "  Source     : " << source_name<< "\n";
    std::cout << "  N_gal      : " << N_gal      << "\n";
    std::cout << "  coord_type : " << coord_type << "\n";
    std::cout << "  boxsize    : " << boxsize    << " cMpc/h\n";
    std::cout << "  redshift   : " << redshift   << "\n";

    // ── Read data ────────────────────────────────────────────────────
    galaxies.reserve(N_gal);
    std::vector<double> buf(4 * N_gal);
    f.read(reinterpret_cast<char*>(buf.data()), N_gal * RECORD_SIZE);

    for (int64_t i = 0; i < N_gal; ++i) {
        Galaxy g;
        g.x      = buf[4*i + 0];
        g.y      = buf[4*i + 1];
        g.z      = buf[4*i + 2];
        g.weight = buf[4*i + 3];
        galaxies.push_back(g);
    }

    // ── Convert coordinates if needed ────────────────────────────────
    if (coord_type == "radec_redshift") {
        std::cout << "  Converting RA/DEC/z → comoving Cartesian...\n";
        convert_radecz(cosmo);
    }

    std::cout << "  Loaded " << galaxies.size() << " galaxies OK\n";
    std::cout << "─────────────────────────────────────────────\n";
}

void GalaxyCatalog::convert_radecz(const BackgroundCosmology &cosmo)
{
    for (auto &g : galaxies) {
        double x_new, y_new, z_new;
        cosmo.radecz_to_xyz(g.x, g.y, g.z, x_new, y_new, z_new);
        g.x = x_new;
        g.y = y_new;
        g.z = z_new;
    }
}

void GalaxyCatalog::load_ascii(const std::string &filename)
{
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("Cannot open: " + filename);

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        Galaxy g;
        g.weight = 1.0;
        ss >> g.x >> g.y >> g.z;
        if (ss.good()) ss >> g.weight;
        galaxies.push_back(g);
    }
    source_name = filename;
    coord_type  = "cartesian_comoving";
    std::cout << "  Loaded " << galaxies.size()
              << " galaxies from ASCII: " << filename << "\n";
}

void GalaxyCatalog::print_summary() const
{
    if (galaxies.empty()) { std::cout << "  Catalog empty\n"; return; }
    double xmin=1e30, xmax=-1e30;
    double ymin=1e30, ymax=-1e30;
    double zmin=1e30, zmax=-1e30;
    for (auto &g : galaxies) {
        xmin=std::min(xmin,g.x); xmax=std::max(xmax,g.x);
        ymin=std::min(ymin,g.y); ymax=std::max(ymax,g.y);
        zmin=std::min(zmin,g.z); zmax=std::max(zmax,g.z);
    }
    std::cout << "\n── Catalog Summary ──────────────────────────\n";
    std::cout << "  N        : " << galaxies.size()  << "\n";
    std::cout << "  x range  : [" << xmin << ", " << xmax << "] Mpc/h\n";
    std::cout << "  y range  : [" << ymin << ", " << ymax << "] Mpc/h\n";
    std::cout << "  z range  : [" << zmin << ", " << zmax << "] Mpc/h\n";
    std::cout << "─────────────────────────────────────────────\n";
}