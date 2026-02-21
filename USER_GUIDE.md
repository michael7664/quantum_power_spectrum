# USER_GUIDE.md — quantum_power_spectrum

This guide documents the `quantum_power_spectrum` project for users who want to run the pipeline and for developers who want to modify or extend it. It covers: what the tool does, how the modules fit together, what files are produced, and where the main “knobs” (parameters/attributes) live in the code.

---

## 1) What this tool does

The repository implements an end-to-end workflow to:

1. Compute a **classical** matter power spectrum \(P_{\rm cl}(k)\) from a baseline cosmology and transfer function.
2. Apply a **QFT-inspired correction** \(\Delta P(k)\) to form \(P_{\rm qu}(k)=P_{\rm cl}(k)+\Delta P(k)\).
3. Transform \(P(k)\) into the real-space two-point correlation function \(\xi(r)\), and show the BAO feature via \(r^2\xi(r)\).
4. Measure \(\xi(r)\) from a tracer catalogue (e.g., TNG-Cluster halos) using the **Landy–Szalay estimator** and compare with the theoretical matter predictions.
5. Produce publication-ready plots from the `.dat` outputs using Python.

Important interpretation note:
- The theory curves produced by this code are for the **matter field**.
- Halo/cluster catalogues are **biased** tracers and are also sparse; their \(\xi(r)\) amplitude and BAO detectability can differ substantially from the matter predictions.

---

## 2) Repository structure (expected)

Typical structure:

```text
quantum_power_spectrum/
├── CMakeLists.txt
├── main.cpp
├── include/
│   ├── BackgroundCosmology.h
│   ├── TransferFunction.h
│   ├── PowerSpectrum.h
│   ├── QuantumQFT.h              (if present; sometimes only .cpp exists)
│   ├── GalaxyCatalog.h
│   ├── LandySzalay.h
│   └── FKPEstimator.h
├── src/
│   ├── BackgroundCosmology.cpp
│   ├── TransferFunction.cpp
│   ├── PowerSpectrum.cpp
│   ├── QuantumQFT.cpp
│   ├── GalaxyCatalog.cpp
│   ├── LandySzalay.cpp
│   └── FKPEstimator.cpp
├── plot_results.py
├── data/
│   └── (catalog files, e.g. TNG-Cluster binary)
└── output/
    ├── *.dat
    └── plots/*.png



3) Build & runtime architecture
3.1 CMake configuration (build-time “knobs”)
File:  CMakeLists.txt 
Key switches/paths you may edit:
	•	C++ standard
	•	 set(CMAKE_CXX_STANDARD 17)  (must be C++17 or newer)
	•	Optimization / debug flags
	•	 CMAKE_CXX_FLAGS_RELEASE  and  CMAKE_CXX_FLAGS_DEBUG 
	•	If debugging numerical issues, prefer Debug builds and add sanitizers if needed.
	•	Quantum++ location
	•	 QPP_DIR  points to  deps/qpp 
	•	Eigen3 include path detection
	•	 find_path(EIGEN3_INCLUDE_DIR ...) 
	•	FFTW3
	•	 find_library(FFTW3_LIB ...)  +  find_path(FFTW3_INCLUDE ...) 
	•	If found, compile definition  HAS_FFTW3  is enabled.
	•	OpenMP
	•	 find_package(OpenMP QUIET) 
	•	If found, links  OpenMP::OpenMP_CXX 
	•	On macOS, there is a Homebrew  libomp  fallback block and a compile definition  HAS_OPENMP=1 .
What to modify in future:
	•	Add options like  -DUSE_FFTW=ON/OFF ,  -DUSE_OPENMP=ON/OFF  to make behavior explicit.
	•	Add  CMAKE_EXPORT_COMPILE_COMMANDS  default ON to help IDE tooling.
3.2 C++ executable flow
File:  main.cpp  (entry point)
Typical responsibilities:
•	Parse CLI arguments (mode selection: theory, data, both).
	•	Load catalogue (if data mode is enabled).
	•	Construct cosmology/transfer function/power spectrum objects.
	•	Call Landy–Szalay measurement if requested.
	•	Write outputs into  output/ .
Developer note:
	•	If you add new outputs, update  plot_results.py  loaders accordingly.


4) Physics modules and editable parameters

4.1 Background cosmology
Files:
	•	 include/BackgroundCosmology.h 
	•	 src/BackgroundCosmology.cpp 
Common attributes (names may vary depending on your implementation):
	•	 H0  or  h 
	•	 Omega_m ,  Omega_b ,  Omega_L ,  Omega_k 
	•	 ns  (spectral index),  As  (amplitude) if you build the primordial spectrum
	•	 sigma8  if you normalize via (\sigma_8)
Where to edit:
•	Look for a parameter struct/class constructor (often in the header).
	•	Look for default values in  BackgroundCosmology.cpp  or where the object is created in  main.cpp .
Future extensions:
	•	Add a config file (YAML/JSON) so users can set cosmology without recompiling.
	•	Add redshift dependence explicitly and store  z  in output headers.

4.2 Transfer function (T(k))
Files:
	•	 include/TransferFunction.h 
	•	 src/TransferFunction.cpp 
Editable features:
	•	Choice of fitting formula and its parameterization.
	•	k-grid definition (range, sampling, spacing: linear/log).
	•	Any baryon wiggle model choices.
Where to edit:
	•	The  TransferFunction  class methods, typically something like  T(k)  or  evaluate(k) .
Future extensions:
	•	Plug in CLASS/CAMB outputs (read transfer function from file).

4.3 Power spectrum (P(k)) and transforms to (\xi(r))
Files:
	•	 include/PowerSpectrum.h 
	•	 src/PowerSpectrum.cpp 
Editable parameters and design points:
	•	k-range:  k_min ,  k_max 
	•	number of k-samples:  N_k 
	•	transform method:
	•	FFTW-based (if  HAS_FFTW3 ), typically faster
	•	direct quadrature/summation fallback, slower but portable
Where to edit:
	•	In  PowerSpectrum.cpp , search for:
	•	k-grid initialization
	•	normalization conventions
	•	Fourier/Hankel transform code path(s)
Future extensions:
	•	Add windowing to reduce ringing in (\xi(r)).
	•	Add output of intermediate integrands for debugging.

4.4 Quantum correction (\Delta P(k))
File:
	•	 src/QuantumQFT.cpp  (and possibly  include/QuantumQFT.h )
Editable features:
	•	The functional form of (\Delta P(k)) (model choice).
	•	Parameters controlling amplitude, scale dependence, IR/UV behavior.
	•	Whether the correction is applied to the primordial spectrum or late-time spectrum.
Where to edit:
	•	Search for code that computes  delta_P ,  dP , or similar arrays.
	•	The output file  power_spectrum_theory_quantum.dat  typically stores columns:
	•	k
	•	(\Delta P(k))
	•	(P_{\rm cl}(k))
	•	(P_{\rm qu}(k))
Future extensions:
	•	Support multiple models with a CLI flag (e.g.  --qft-model modelA|modelB ).
	•	Add parameter priors and scanning tools.


    5) Data ingestion: GalaxyCatalog
Files:
	•	 include/GalaxyCatalog.h 
	•	 src/GalaxyCatalog.cpp 
Typical attributes you will find (based on your earlier code style):
	•	 galaxies : a vector of objects with  (x, y, z, w) 
	•	metadata:
	•	 source_name 
	•	 coord_type  (e.g.  "cartesian_comoving" )
	•	 boxsize 
	•	 redshift 
	•	 is_periodic 
Editable features:
	•	Catalogue reader: binary format decoding and validation.
	•	Unit conventions (Mpc/h vs Mpc).
	•	Periodic boundary handling (if supported by your pair counter; see LandySzalay section).
Where to edit:
•	In  GalaxyCatalog.cpp , locate the loader for your catalogue path.
	•	If you add extra fields (mass, ID), keep them optional and ensure the reader remains backward-compatible.
Future extensions:
	•	Support ASCII/CSV and HDF5 inputs.
	•	Add selection cuts (mass threshold, region mask) as CLI options.


    6) Landy–Szalay correlation function measurement
Files:
	•	 include/LandySzalay.h 
	•	 src/LandySzalay.cpp 
6.1 Core estimator
The Landy–Szalay estimator is:
\xi_{\rm LS}(r)=\frac{DD(r)-2DR(r)+RR(r)}{RR(r)}.
Implementation note (critical):
	•	If your code uses normalized pair fractions (recommended), you must normalize each of DD, DR, RR by the corresponding total number of possible pairs:
	•	(n_{DD}=N_D(N_D-1)/2)
	•	(n_{DR}=N_DN_R)
	•	(n_{RR}=N_R(N_R-1)/2)

6.2 Editable parameters (“knobs”) in LandySzalay

(A) Radial binning
	•	Location:  LandySzalay::LandySzalay(double r_min, double r_max, int N_bins) 
	•	Parameters:
	•	 r_min  : minimum separation included
	•	 r_max  : maximum separation included
	•	 N_bins : number of logarithmic bins
Where to edit:
	•	Wherever the  LandySzalay  object is constructed (often in  main.cpp ).
	•	The bin mapping is typically in  bin_index(r) .

(B) Random catalogue size & seed
	•	Location:  LandySzalay::make_randoms(const GalaxyCatalog&, int N_randoms, int seed) 
	•	Parameters:
	•	 N_randoms : if negative, the code uses a default multiple of  N_data  (commonly 5×)
	•	 seed : RNG seed for reproducibility
Where to edit:
	•	The call site (again typically  main.cpp ) or the default rule in  make_randoms .
Recommended practice:
	•	For stable (\xi(r)), use  N_randoms >= 20 × N_data  when feasible.

(C) Geometry / selection function for randoms
	•	Current simple approach (common in prototypes): uniform randoms inside the data bounding box.
	•	Editable to:
	•	periodic-box randoms
	•	survey mask / radial selection
Where to edit:
	•	 make_randoms : change how  (x,y,z)  are sampled.

(D) Pair counting
	•	Location:  LandySzalay::count_pairs(A, B, autocorr) 
	•	Editable features:
	•	periodic boundary conditions (if not implemented, adding minimal-image convention is a common next step)
	•	performance: OpenMP scheduling, chunk size, SIMD, cell lists, kd-trees
Where to edit:
	•	The loops inside  count_pairs .
	•	OpenMP blocks are guarded by  #ifdef HAS_OPENMP .

(E) Error model
	•	Common quick estimate: Poisson  err ~ (1+xi)/sqrt(DD) 
	•	Future improvements:
	•	jackknife resampling
	•	bootstrap resampling
	•	mocks / covariance estimation
Where to edit:
	•	In  compute() , the error assignment line.


6.3 Output file format:  xi_data_classical.dat 
File written by  LandySzalay::write(...)  typically follows:

# r[Mpc/h]  xi_LS(r)  error  DD  DR  RR
r  xi  err  DD  DR  RR
...


Editable features:
	•	Add more columns (e.g. normalized DD/DR/RR, or  n_DD  etc.) for debugging.
	•	Add a header recording  ND ,  NR ,  r_min ,  r_max ,  N_bins ,  seed .
Where to edit:
	•	 LandySzalay::write() .


7) Plotting: plot_results.py
File:  plot_results.py 

7.1 What the script expects to find
It loads  .dat  files from the C++ pipeline and writes PNG figures. Typical expected inputs:
	•	 output/power_spectrum_theory.dat 
	•	 output/power_spectrum_theory_quantum.dat 
	•	 output/power_spectrum_data.dat 
	•	 output/xi_theory_classical.dat 
	•	 output/xi_theory_quantum.dat 
	•	 output/xi_data_classical.dat 
	•	(optional)  output/xi_data_quantum.dat ,  output/power_spectrum_data_quantum.dat 
If you rename output files or move them to subdirectories, update the loader paths in  plot_results.py .

7.2 Script-level editable parameters

(A) Matplotlib style
	•	Location:  plt.rcParams.update({...}) 
	•	Edit:
	•	font family, sizes, DPI, grid style, etc.
(B) Color palette
	•	Location:  COLORS = {...} 
	•	Edit to match your preferred publication style.
(C) Smoothing / ringing control
	•	Location:
	•	 load_xi_theory_qu(...)  and often  load_xi_data_qu(...) 
	•	The script currently applies a small moving average via  uniform_filter1d  to suppress Fourier ringing.
To change behavior:
	•	Disable smoothing: return  xi  directly.
	•	Change smoothing strength: edit the window size ( size ).
	•	Replace with Savitzky–Golay filter (preserves peaks better).
(D) Column conventions
	•	Location: loader functions such as:
	•	 load_pk_quantum(...) 
	•	 load_xi_data_cl(...) 
If your C++ output format changes (columns added/reordered), update:
	•	 usecols=(...) 
	•	the mapping logic to  k ,  dP ,  P_cl , etc.
(E) Plot ranges
	•	Locations:  ax.set_xlim(...) ,  ax.set_ylim(...) , etc.
	•	Edit for:
	•	k-range shown in (P(k))
	•	r-range shown in (\xi(r))
	•	BAO window limits for  r^2\xi(r)  panel
(F) Legend labels
	•	If you want to hide data legend entries (e.g., for a “theory-only legend”), remove  label=...  for the TNG data series calls.


    7.3 CLI flags and output folder
The script uses an argument for where to write plots (named  --output-dir  in the current version). file:174
Where to edit:
	•	The  argparse  section near the bottom of the script.
	•	The calls to  os.makedirs(...)  and  savefig(...) .
Future extension:
	•	Add  --input-dir  to read  .dat  from a non-default location.
	•	Add  --no-smooth  to disable smoothing from the CLI.
	•	Add  --format pdf|png|svg  to export vector graphics.


8) Output files (what they mean)

8.1 Theory power spectrum files
 output/power_spectrum_theory.dat 
	•	Intended columns:
	•	 k h/Mpc 
	•	 P_classical(k) (Mpc/h)^3 
 output/power_spectrum_theory_quantum.dat 
	•	Intended columns:
	•	 k 
	•	 delta_P 
	•	 P_classical 
	•	 P_quantum  (or a column from which  P_quantum  can be reconstructed)
Developer note:
	•	If you store  delta_P  and  P_classical , you can always reconstruct  P_quantum  robustly.


8.2 Theory correlation functions
 output/xi_theory_classical.dat 
	•	Intended columns:
	•	 r Mpc/h 
	•	 xi_classical(r) 
 output/xi_theory_quantum.dat 
	•	Intended columns:
	•	 r 
	•	 xi_quantum(r) 

8.3 Data correlation function (Landy–Szalay)
 output/xi_data_classical.dat 
	•	Intended columns:
	•	 r 
	•	 xi_LS(r) 
	•	 error 
	•	 DD ,  DR ,  RR 


9) Extension points (recommended roadmap)

9.1 Improve random catalogue realism
Current randoms may be uniform in a bounding box. For more realistic survey-like conditions:
	•	implement a selection function (\bar{n}(\vec{x}))
	•	apply survey mask and radial distribution
	•	incorporate periodic boundary conditions if your simulation is periodic
Where:
	•	 LandySzalay::make_randoms() 

9.2 Add periodic boundary conditions (PBC)
If your simulation volume is periodic and you want correct large-scale pair distances:
	•	use minimum-image convention for  dx,dy,dz 
Where:
	•	 LandySzalay::count_pairs() 

9.3 Replace O(N²) pair counting with an accelerator
For large catalogues:
	•	grid / cell lists
	•	kd-trees / ball trees
	•	dual-tree pair counting
Where:
	•	new class (e.g.  PairCounter ) called by  LandySzalay .

9.4 Add covariance estimation
For robust BAO inference:
	•	jackknife regions (split box into subvolumes)
	•	bootstrap resampling
	•	mock catalogues
Where:
	•	add an analysis module + store covariance output (e.g.,  cov_xi.dat )

9.5 Model halo bias explicitly (optional)
If you want direct comparison of matter theory to halo measurements:
	•	estimate bias (b) from large-scale ratio
	•	apply (b^2) scaling to theory curves (in plot script or analysis step)
	•	optionally allow scale-dependent bias
Where:
	•	analysis layer (Python) is usually best.


0) Troubleshooting and validation checklist

10.1 Sanity checks for (\xi(r)) output
	•	DD/DR/RR should increase with r-bin volume at large r (roughly).
	•	(\xi(r)) should not be a nearly constant ~1 value across all r (that usually indicates a normalization bug).
	•	At large r, many bins may have DD=0; ensure you handle them (masking or setting xi=0).

10.2 Plot mismatches (theory vs TNG clusters)
	•	A mismatch in amplitude between matter theory and cluster measurements is expected (bias).
	•	BAO non-detection is expected for small, sparse samples; don’t over-interpret.

10.3 Performance issues
	•	Enable OpenMP if possible (see CMake output).
	•	Increase  N_randoms  only as far as needed for stability; it increases DR and RR cost.



11) Developer notes: documenting parameters
If you plan to share this tool publicly, consider adding:
	•	a single  config.yaml  with:
	•	cosmology parameters
	•	k-grid and r-grid settings
	•	Landy–Szalay settings ( r_min ,  r_max ,  N_bins ,  N_randoms ,  seed )
	•	a  --config  flag in  main.cpp 
    •	output headers that echo the runtime configuration into each  .dat  file
This will make runs reproducible without relying on commit history.


12) Contact / contributions
If you extend the code (new QFT models, improved estimators, new catalogue formats), please:
	•	document the new parameters and defaults,
	•	keep output formats stable or versioned,
	•	add short examples to the README showing new flags.