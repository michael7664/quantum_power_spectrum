cmake_minimum_required(VERSION 3.16)
project(quantum_power_spectrum VERSION 1.0 LANGUAGES CXX)

# ── C++ Standard ──────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Build type default ────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler   : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# ── Compiler flags ────────────────────────────────────────────────────
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -Wall -Wextra -Wpedantic")

# ── Quantum++ (header-only) ───────────────────────────────────────────
set(QPP_DIR "${CMAKE_SOURCE_DIR}/deps/qpp")

if(NOT EXISTS "${QPP_DIR}/include/qpp/qpp.hpp")
    message(FATAL_ERROR
        "Quantum++ not found at ${QPP_DIR}/include/qpp/qpp.hpp\n"
        "Run:  git clone --depth=1 https://github.com/softwareQinc/qpp.git deps/qpp")
endif()

message(STATUS "Quantum++  : ${QPP_DIR}")

# ── Quantum++ requires Eigen3 ─────────────────────────────────────────
# Search all common locations inside qpp
find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Dense
    HINTS
        "${QPP_DIR}/include"
        "${QPP_DIR}/include/eigen3"
        "${QPP_DIR}/eigen3"
        "${QPP_DIR}/deps/eigen3"
        /opt/homebrew/include/eigen3
        /usr/local/include/eigen3
    NO_DEFAULT_PATH
)

if(NOT EIGEN3_INCLUDE_DIR)
    # Fall back to system Eigen via find_package
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIRS}")
        message(STATUS "Eigen3     : system (${EIGEN3_INCLUDE_DIR})")
    else()
        message(FATAL_ERROR
            "Eigen3 not found.\n"
            "Install via:  brew install eigen\n"
            "Or check:     find deps/qpp -name Eigen -type d")
    endif()
else()
    message(STATUS "Eigen3     : ${EIGEN3_INCLUDE_DIR}")
endif()

# ── Find OpenMP (optional but recommended for pair counting) ──────────
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP     : found (${OpenMP_CXX_VERSION})")
else()
    message(STATUS "OpenMP     : not found — single-threaded pair counting")
endif()

# ── Find FFTW3 (optional — used by PowerSpectrum.cpp) ─────────────────
find_library(FFTW3_LIB fftw3
    HINTS /usr/local/lib /opt/homebrew/lib /usr/lib)
find_path(FFTW3_INCLUDE fftw3.h
    HINTS /usr/local/include /opt/homebrew/include /usr/include)

if(FFTW3_LIB AND FFTW3_INCLUDE)
    message(STATUS "FFTW3      : found (${FFTW3_LIB})")
    set(HAS_FFTW3 TRUE)
else()
    message(STATUS "FFTW3      : not found — using direct summation fallback")
    set(HAS_FFTW3 FALSE)
endif()

# ── Source files ──────────────────────────────────────────────────────
set(SOURCES
    src/BackgroundCosmology.cpp
    src/TransferFunction.cpp
    src/PowerSpectrum.cpp
    src/QuantumQFT.cpp
    src/GalaxyCatalog.cpp
    src/LandySzalay.cpp
    src/FKPEstimator.cpp
)

# ── Executable ────────────────────────────────────────────────────────
add_executable(quantum_ps main.cpp ${SOURCES})

# ── Include directories ───────────────────────────────────────────────
target_include_directories(quantum_ps PRIVATE
    ${CMAKE_SOURCE_DIR}/include          # project headers
    ${QPP_DIR}/include                   # Quantum++ top-level
    ${EIGEN3_INCLUDE_DIR}                # Eigen3
)

# ── Compile definitions ───────────────────────────────────────────────
if(HAS_FFTW3)
    target_compile_definitions(quantum_ps PRIVATE HAS_FFTW3=1)
    target_include_directories(quantum_ps PRIVATE ${FFTW3_INCLUDE})
endif()

# ── Link libraries ────────────────────────────────────────────────────
target_link_libraries(quantum_ps PRIVATE
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${HAS_FFTW3}>:${FFTW3_LIB}>
    m      # math library
)

# ── macOS specific: suppress OpenMP warning from Apple Clang ──────────
if(APPLE)
    if(NOT OpenMP_CXX_FOUND)
        # Try libomp from Homebrew explicitly
        find_library(OMP_LIB omp
            HINTS /opt/homebrew/opt/libomp/lib
                  /usr/local/opt/libomp/lib)
        find_path(OMP_INCLUDE omp.h
            HINTS /opt/homebrew/opt/libomp/include
                  /usr/local/opt/libomp/include)
        if(OMP_LIB AND OMP_INCLUDE)
            message(STATUS "OpenMP     : found via Homebrew libomp")
            target_include_directories(quantum_ps PRIVATE ${OMP_INCLUDE})
            target_link_libraries(quantum_ps PRIVATE ${OMP_LIB})
            target_compile_options(quantum_ps PRIVATE -Xpreprocessor -fopenmp)
            target_compile_definitions(quantum_ps PRIVATE HAS_OPENMP=1)
        endif()
    endif()
    # Suppress LLVM/Clang deprecation warnings on macOS
    target_compile_options(quantum_ps PRIVATE -Wno-deprecated-declarations)
endif()

# ── Output directories ────────────────────────────────────────────────
set_target_properties(quantum_ps PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build"
)

# Make sure output/ directory exists at build time
add_custom_command(TARGET quantum_ps POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_SOURCE_DIR}/output"
    COMMENT "Ensuring output/ directory exists"
)

# ── Print final summary ───────────────────────────────────────────────
message(STATUS "")
message(STATUS "══════════════════════════════════════════")
message(STATUS "  quantum_power_spectrum build config")
message(STATUS "══════════════════════════════════════════")
message(STATUS "  Sources    : ${SOURCES}")
message(STATUS "  OpenMP     : ${OpenMP_CXX_FOUND}")
message(STATUS "  FFTW3      : ${HAS_FFTW3}")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "══════════════════════════════════════════")
message(STATUS "")
