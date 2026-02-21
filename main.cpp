#include <iostream>
#include <string>
#include <stdexcept>
#include <numeric>
#include <cmath>

#include "BackgroundCosmology.h"
#include "TransferFunction.h"
#include "PowerSpectrum.h"
#include "QuantumQFT.h"
#include "GalaxyCatalog.h"
#include "LandySzalay.h"
#include "FKPEstimator.h"

// ── Estimate median halo mass from catalog source name ────────────────
static double estimate_median_mass(const GalaxyCatalog &cat) {
    if (cat.source_name.find("TNG-Cluster") != std::string::npos)
        return 7.0e14;   // Msun — median of our 352-cluster sample
    if (cat.source_name.find("TNG300") != std::string::npos
     || cat.source_name.find("TNG100") != std::string::npos)
        return 1.0e12;
    return 3.0e12;
}

// ── Theory mode ───────────────────────────────────────────────────────
static void run_theory(const BackgroundCosmology &cosmo,
                       const TransferFunction    &tf,
                       double b_halo = 1.0)
{
    std::cout << "\n════════════════════════════════════════════\n";
    std::cout << "  MODE: THEORY  (b_halo = " << b_halo << ")\n";
    std::cout << "════════════════════════════════════════════\n";

    cosmo.print();

    PowerSpectrum ps(cosmo, tf);
    QuantumQFT    qft(cosmo);

    // Matter P(k) (b=1) and biased P(k) (b=b_halo)
    auto pk_matter = ps.compute(1e-3, 10.0, 4000, 1.0);
    auto pk_biased = ps.compute(1e-3, 10.0, 4000, b_halo);
    PowerSpectrum::write(pk_matter, "output/power_spectrum_theory.dat");
    PowerSpectrum::write(pk_biased, "output/power_spectrum_theory_biased.dat");

    // Quantum-corrected P(k) from biased spectrum
    std::vector<std::pair<double,double>> pk_pairs;
    pk_pairs.reserve(pk_biased.size());
    for (auto &p : pk_biased)
        pk_pairs.push_back({p.k, p.Pk});
    auto qft_corr = qft.apply(pk_pairs);
    QuantumQFT::write(qft_corr, "output/power_spectrum_theory_quantum.dat");

    // Theory ξ(r) — linear r-grid, 1–250 Mpc/h, 100 points
    auto xi_cl = ps.xi_theory(1.0, 250.0, 1000, b_halo);
    PowerSpectrum::write_xi(xi_cl, "output/xi_theory_classical.dat");

    auto xi_qu = qft.xi_quantum(qft_corr, 1.0, 250.0, 1000);
    QuantumQFT::write_xi(xi_qu, "output/xi_theory_quantum.dat");

    std::cout << "\n── Theory outputs written to output/ ────────\n";
    std::cout << "  b_halo = " << b_halo
              << "  →  b² = " << b_halo*b_halo << "\n";
}

// ── Data mode ─────────────────────────────────────────────────────────
static void run_data(const BackgroundCosmology &cosmo,
                     const TransferFunction    &tf,
                     const std::string         &bin_file)
{
    std::cout << "\n════════════════════════════════════════════\n";
    std::cout << "  MODE: DATA  →  " << bin_file << "\n";
    std::cout << "════════════════════════════════════════════\n";

    GalaxyCatalog data;
    data.load_bin(bin_file, cosmo);
    data.print_summary();

    // Estimate halo bias from median mass (Tinker et al. 2010)
    double M_median = estimate_median_mass(data);
    double b_halo   = PowerSpectrum::bias_from_mass(
                          M_median, data.redshift);
    std::cout << "\n  Median M200c : " << M_median << " Msun\n";
    std::cout << "  Linear bias  : b = " << b_halo << "\n";
    std::cout << "  P(k) boost   : b² = "
              << b_halo*b_halo << "\n";

    // 50× randoms for reliable Landy-Szalay
    // (minimum ~17,600 for 352 clusters)
    int N_randoms = 50 * (int)data.size();
    std::cout << "  Randoms      : " << N_randoms << "\n";

    // ── Landy-Szalay ξ(r) ────────────────────────────────────────
    std::cout << "\n── Landy-Szalay pair counting ───────────────\n";
    // r range: 5–400 Mpc/h  (covers full box separation range)
    // 25 log-spaced bins
    LandySzalay ls(5.0, 400.0, 25);
    auto randoms = ls.make_randoms(data, N_randoms);
    auto xi_ls   = ls.compute(data, randoms);
    LandySzalay::write(xi_ls, "output/xi_data_classical.dat");

    // ── FKP P(k) via FFT grid ─────────────────────────────────────
    std::cout << "\n── FKP P(k) estimator ───────────────────────\n";
    // Ng=64 grid, k range 0.005–0.5 h/Mpc, 30 bins
    FKPEstimator fkp(1e4, 0.005, 0.5, 30, 64);
    auto pk_data = fkp.compute(data, randoms, cosmo);
    FKPEstimator::write(pk_data, "output/power_spectrum_data.dat");

    // ── Quantum corrections on data P(k) ─────────────────────────
    QuantumQFT qft(cosmo);
    std::vector<std::pair<double,double>> pk_pairs;
    for (auto &p : pk_data)
        pk_pairs.push_back({p.k, p.Pk_shot_subtracted});
    auto qft_corr = qft.apply(pk_pairs);
    QuantumQFT::write(qft_corr, "output/power_spectrum_data_quantum.dat");

    // Quantum ξ(r) from data P(k) — linear r-grid 5–400 Mpc/h
    auto xi_qu = qft.xi_quantum(qft_corr, 5.0, 400.0, 25);
    QuantumQFT::write_xi(xi_qu, "output/xi_data_quantum.dat");

    std::cout << "\n── Data outputs written to output/ ──────────\n";
    std::cout << "  bias used: b = " << b_halo << "\n";
}

// ── Main ──────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    std::cout << "\n╔══════════════════════════════════════════╗\n";
    std::cout << "║   Quantum Power Spectrum  v1.2           ║\n";
    std::cout << "╚══════════════════════════════════════════╝\n";

    if (argc < 2) {
        std::cout << "\nUsage:\n"
                  << "  ./quantum_ps theory\n"
                  << "  ./quantum_ps tng   <catalog.bin>\n"
                  << "  ./quantum_ps both  <catalog.bin>\n\n";
        return 1;
    }

    std::string mode     = argv[1];
    std::string bin_file = (argc >= 3) ? argv[2] : "";

    try {
        BackgroundCosmology cosmo;
        TransferFunction    tf(cosmo);

        if (mode == "theory") {
            double b = PowerSpectrum::bias_from_mass(7.0e14, 0.0);
            run_theory(cosmo, tf, b);

        } else if (mode == "tng"  ||
                   mode == "sdss" ||
                   mode == "data") {
            if (bin_file.empty())
                throw std::runtime_error(
                    "Usage: ./quantum_ps " + mode + " <catalog.bin>");
            run_data(cosmo, tf, bin_file);

        } else if (mode == "both") {
            if (bin_file.empty())
                throw std::runtime_error(
                    "Usage: ./quantum_ps both <catalog.bin>");
            // Load catalog once to determine bias
            GalaxyCatalog tmp;
            tmp.load_bin(bin_file, cosmo);
            double M_med  = estimate_median_mass(tmp);
            double b_halo = PowerSpectrum::bias_from_mass(
                                M_med, tmp.redshift);
            run_theory(cosmo, tf, b_halo);
            run_data  (cosmo, tf, bin_file);

        } else {
            throw std::runtime_error(
                "Unknown mode: " + mode +
                "\nValid: theory | tng | sdss | data | both");
        }

    } catch (const std::exception &e) {
        std::cerr << "\n[ERROR] " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n╔══════════════════════════════════════════╗\n";
    std::cout << "║   Done. Results in output/               ║\n";
    std::cout << "╚══════════════════════════════════════════╝\n\n";
    return 0;
}
