#!/usr/bin/env python3
"""
fetch_sdss.py
─────────────
Download SDSS DR17 galaxy catalog via astroquery and write a
universal .bin for the quantum_power_spectrum C++ project.

Coordinates are stored as (RA, DEC, redshift) — the C++ project
reads coord_type="radec_redshift" and converts via BackgroundCosmology.

Requires:
  pip install astroquery astropy numpy

Usage:
  python fetch_sdss.py [options]

Options:
  --out     Output .bin path          (default: auto-named)
  --z-min   Min redshift              (default: 0.02)
  --z-max   Max redshift              (default: 0.15)
  --N-max   Max galaxies              (default: 50000)
  --inspect Inspect existing .bin and exit
"""

import argparse, os
import numpy as np
from catalog_io import write_catalog, inspect_catalog


def fetch_sdss_galaxies(z_min, z_max, N_max):
    try:
        from astroquery.sdss import SDSS
    except ImportError:
        print("ERROR: astroquery not found.  pip install astroquery astropy")
        raise

    query = f"""
    SELECT TOP {N_max}
        p.ra,
        p.dec,
        s.z            AS redshift,
        s.zErr         AS z_err,
        1.0            AS weight
    FROM PhotoObj AS p
    JOIN SpecObj  AS s ON s.bestobjid = p.objid
    WHERE s.class      = 'GALAXY'
      AND s.zWarning   = 0
      AND s.z          BETWEEN {z_min} AND {z_max}
      AND p.clean      = 1
      AND p.mode       = 1
    ORDER BY NEWID()
    """

    print(f"Querying SDSS DR17 (z=[{z_min},{z_max}], N_max={N_max})...")
    result = SDSS.query_sql(query, data_release=17, timeout=300)

    if result is None or len(result) == 0:
        raise RuntimeError("SDSS query returned no results")

    ra  = np.array(result["ra"],       dtype=np.float64)
    dec = np.array(result["dec"],      dtype=np.float64)
    z   = np.array(result["redshift"], dtype=np.float64)
    w   = np.ones(len(ra),             dtype=np.float64)

    print(f"  Fetched {len(ra):,} galaxies")
    print(f"  RA  range : [{ra.min():.2f},  {ra.max():.2f}] deg")
    print(f"  DEC range : [{dec.min():.2f}, {dec.max():.2f}] deg")
    print(f"  z   range : [{z.min():.4f},  {z.max():.4f}]")

    positions = np.column_stack([ra, dec, z])
    return positions, w, float(np.median(z))


def main():
    p = argparse.ArgumentParser(
        description="Fetch SDSS DR17 → universal .bin",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    p.add_argument("--out",     default=None)
    p.add_argument("--z-min",   type=float, default=0.02)
    p.add_argument("--z-max",   type=float, default=0.15)
    p.add_argument("--N-max",   type=int,   default=50000)
    p.add_argument("--inspect", default=None)
    args = p.parse_args()

    if args.inspect:
        inspect_catalog(args.inspect)
        return

    if args.out is None:
        args.out = (f"../quantum_power_spectrum/data/"
                    f"SDSS_DR17_z{args.z_min:.2f}-{args.z_max:.2f}.bin")

    positions, weights, z_eff = fetch_sdss_galaxies(
        args.z_min, args.z_max, args.N_max)

    write_catalog(
        filename    = args.out,
        positions   = positions,
        weights     = weights,
        boxsize     = 0.0,           # not a periodic box
        redshift    = z_eff,
        is_periodic = False,         # survey — needs mask + randoms
        has_weights = False,
        source_name = f"SDSS_DR17_z{args.z_min:.2f}-{args.z_max:.2f}",
        coord_type  = "radec_redshift",  # C++ converts RA,DEC,z → x,y,z
    )
    inspect_catalog(args.out)


if __name__ == "__main__":
    main()
