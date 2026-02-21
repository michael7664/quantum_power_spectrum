#!/usr/bin/env python3
"""
convert_custom.py
─────────────────
Convert any user-supplied catalog (ASCII, CSV, FITS) to the
universal .bin format for the quantum_power_spectrum C++ project.

Supports:
  .txt / .dat  — whitespace-separated, first line may be header
  .csv         — comma-separated
  .fits / .fits.gz — FITS binary table

Usage examples:

  # Cartesian comoving (x y z in cMpc/h), periodic simulation box:
  python convert_custom.py \
      --in my_sim.txt --cols x y z \
      --coord-type cartesian_comoving \
      --boxsize 500 --redshift 0.0 --periodic \
      --source-name "MyNbody_z0" \
      --out data/my_sim.bin

  # Cartesian with weights:
  python convert_custom.py \
      --in my_sim.csv --cols x y z weight \
      --coord-type cartesian_comoving \
      --boxsize 200 --redshift 1.0 --periodic \
      --out data/my_sim_weighted.bin

  # Angular survey (RA DEC redshift):
  python convert_custom.py \
      --in survey.fits --cols RA DEC Z_SPEC \
      --coord-type radec_redshift \
      --boxsize 0 --redshift 0.3 \
      --source-name "MySurvey" \
      --out data/my_survey.bin

  # Inspect an existing .bin:
  python convert_custom.py --inspect data/my_sim.bin
"""

import argparse, os, sys
import numpy as np
from catalog_io import write_catalog, inspect_catalog


def load_fits(filepath, col_names):
    try:
        from astropy.table import Table
    except ImportError:
        print("ERROR: astropy required for FITS.  pip install astropy")
        sys.exit(1)

    t = Table.read(filepath)
    cols = []
    for c in col_names:
        if c not in t.colnames:
            raise KeyError(f"Column '{c}' not found in FITS. "
                           f"Available: {t.colnames}")
        cols.append(np.array(t[c], dtype=np.float64))

    positions = np.column_stack(cols[:3])
    weights   = cols[3] if len(cols) > 3 else np.ones(len(positions))
    return positions, weights


def load_text(filepath, col_names):
    """
    Load plain text / CSV.
    Tries to auto-detect column names from first comment line.
    Falls back to integer column indices.
    """
    delimiter = "," if filepath.lower().endswith(".csv") else None

    # Try to read header
    col_index = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                header_cols = line.lstrip("#").split(delimiter or None)
                col_index = {c.strip(): i for i, c in enumerate(header_cols)}
                break

    data = np.loadtxt(filepath, comments="#", delimiter=delimiter)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    def get_col(name):
        if name in col_index:
            return data[:, col_index[name]]
        try:
            return data[:, int(name)]
        except (ValueError, IndexError):
            raise KeyError(f"Column '{name}' not found. "
                           f"Header columns: {list(col_index.keys())}")

    pos_cols  = [get_col(c) for c in col_names[:3]]
    positions = np.column_stack(pos_cols)
    weights   = get_col(col_names[3]) if len(col_names) > 3 \
                else np.ones(len(positions))
    return positions, weights


def main():
    p = argparse.ArgumentParser(
        description="Convert any catalog → universal .bin",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    p.add_argument("--in",          dest="infile",  required=True,
                   help="Input file (.txt, .csv, .dat, .fits, .fits.gz)")
    p.add_argument("--out",         required=True,
                   help="Output .bin file path")
    p.add_argument("--cols",        nargs="+", required=True,
                   help="Column names or indices: x y z [weight]")
    p.add_argument("--coord-type",  default="cartesian_comoving",
                   choices=["cartesian_comoving", "radec_redshift"],
                   help="Coordinate type")
    p.add_argument("--boxsize",     type=float, default=0.0,
                   help="Box size in cMpc/h (0 for surveys)")
    p.add_argument("--redshift",    type=float, default=0.0,
                   help="Effective redshift of the sample")
    p.add_argument("--periodic",    action="store_true",
                   help="Flag if simulation has periodic boundaries")
    p.add_argument("--source-name", default=None,
                   help="Label for the catalog (default: input filename)")
    p.add_argument("--N-max",       type=int, default=None,
                   help="Maximum number of objects to keep")
    p.add_argument("--seed",        type=int, default=42)
    p.add_argument("--inspect",     default=None,
                   help="Inspect existing .bin file and exit")
    args = p.parse_args()

    if args.inspect:
        inspect_catalog(args.inspect)
        return

    np.random.seed(args.seed)

    # ── Load input file ───────────────────────────────────────────────
    ext = args.infile.lower()
    print(f"Loading: {args.infile}")

    if ext.endswith(".fits") or ext.endswith(".fits.gz"):
        positions, weights = load_fits(args.infile, args.cols)
    else:
        positions, weights = load_text(args.infile, args.cols)

    print(f"  Loaded {len(positions):,} objects")

    # ── Subsample ─────────────────────────────────────────────────────
    if args.N_max and len(positions) > args.N_max:
        idx       = np.random.choice(len(positions), args.N_max, replace=False)
        positions = positions[idx]
        weights   = weights[idx]
        print(f"  Subsampled to {args.N_max:,}")

    # ── Quick stats ───────────────────────────────────────────────────
    for i, ax in enumerate(["col0", "col1", "col2"]):
        lo, hi = positions[:, i].min(), positions[:, i].max()
        print(f"  {args.cols[i]}: [{lo:.4g}, {hi:.4g}]")

    # ── Source name ───────────────────────────────────────────────────
    src = args.source_name or os.path.splitext(
              os.path.basename(args.infile))[0]

    # ── Write .bin ────────────────────────────────────────────────────
    write_catalog(
        filename    = args.out,
        positions   = positions,
        weights     = weights,
        boxsize     = args.boxsize,
        redshift    = args.redshift,
        is_periodic = args.periodic,
        has_weights = len(args.cols) > 3,
        source_name = src,
        coord_type  = args.coord_type,
    )
    inspect_catalog(args.out)


if __name__ == "__main__":
    main()
