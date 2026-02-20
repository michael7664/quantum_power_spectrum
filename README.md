# quantum_power_spectrum

End‑to‑end pipeline to compute quantum field theory corrections to the cosmological power spectrum and correlation function.

---

## 1. Requirements

- CMake ≥ 3.16
- C++17 compiler
- Python 3.9+ (`numpy`, `matplotlib`, `scipy`)
- Optional: FFTW3, OpenMP

---

## 2. Build

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j$(nproc)
cd ..
./build# Quantum-Power-Spectrum
