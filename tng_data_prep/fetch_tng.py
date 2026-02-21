#!/usr/bin/env python3
"""
fetch_tng.py
────────────
Download TNG50 / TNG100 / TNG300 subhalo catalog and write a
universal .bin for the quantum_power_spectrum C++ project.

These are PERIODIC uniform boxes (unlike TNG-Cluster).

Usage:
  python fetch_tng.py --api-key $TNG_API_KEY [options]

Options:
  --sim       TNG50-1 | TNG100-1 | TNG300-1   (default: TNG100-1)
  --snap      Snapshot 0-99                   (default: 99 = z=0)
  --mass-min  Min stellar mass [Msun]         (default: 1e9)
  --N-max     Max galaxies                    (default: 100000)
  --method    bulk (fast) | paged (slow)      (default: bulk)
  --out       Output .bin path                (auto-named if omitted)
  --inspect   Inspect existing .bin and exit
  --dry-run   Print plan, no download
"""

import argparse, os, sys, io, time
import numpy as np
import requests

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

from catalog_io import write_catalog, inspect_catalog

BASE_URL  = "https://www.tng-project.org/api"
TNG_H     = 0.6774
MASS_UNIT = 1e10 / TNG_H    # code mass → Msun
POS_UNIT  = 1.0  / 1000.0   # ckpc/h   → cMpc/h

TNG_BOXES = {
    "TNG50-1":   35.0,  "TNG50-2":   35.0,
    "TNG50-3":   35.0,  "TNG50-4":   35.0,
    "TNG100-1":  75.0,  "TNG100-2":  75.0,  "TNG100-3":  75.0,
    "TNG300-1": 205.0,  "TNG300-2": 205.0,  "TNG300-3": 205.0,
}


def api_get(url, api_key, params=None, binary=False, retries=3):
    headers = {"api-key": api_key}
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers,
                             params=params, timeout=180)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 15))
                print(f"  Rate limited — waiting {wait}s...")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.content if binary else r.json()
        except requests.exceptions.Timeout:
            print(f"  Timeout attempt {attempt+1}/{retries}...")
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Failed: {url}")


def get_snapshot_redshift(api_key, sim, snap):
    data = api_get(f"{BASE_URL}/{sim}/snapshots/{snap}/", api_key)
    return data["redshift"]


# ── Bulk method: one HDF5 download ───────────────────────────────────
def fetch_bulk(api_key, sim, snap, mass_min_msun, N_max, dry_run=False):
    if not HAS_H5PY:
        print("ERROR: h5py required.  pip install h5py")
        sys.exit(1)

    url    = f"{BASE_URL}/{sim}/snapshots/{snap}/groupcat/"
    params = {"Subhalo": "SubhaloPos,SubhaloMassType"}

    if dry_run:
        print(f"  [DRY RUN] GET {url}  params={params}")
        return None

    print(f"  Downloading subhalo catalog for {sim} snap={snap}...")
    raw  = api_get(url, api_key, params=params, binary=True)
    print(f"  Downloaded {len(raw)/1e6:.1f} MB")

    with h5py.File(io.BytesIO(raw), "r") as f:
        pos      = f["Subhalo/SubhaloPos"][:]           # ckpc/h
        mass_all = f["Subhalo/SubhaloMassType"][:]      # 1e10 Msun/h

    m_star = mass_all[:, 4] * MASS_UNIT                 # stellar → Msun
    mask   = m_star >= mass_min_msun
    pos    = pos[mask]
    m_star = m_star[mask]
    print(f"  Subhalos after M* > {mass_min_msun:.1e} Msun : {len(pos)}")

    # Subsample if needed
    if len(pos) > N_max:
        idx    = np.random.choice(len(pos), N_max, replace=False)
        pos    = pos[idx]
        m_star = m_star[idx]
        print(f"  Subsampled to {N_max}")

    pos_mpc = pos * POS_UNIT
    print(f"  Position range x: [{pos_mpc[:,0].min():.1f},"
          f" {pos_mpc[:,0].max():.1f}] cMpc/h")
    return pos_mpc


# ── Paged method: REST API pagination (no h5py needed) ───────────────
def fetch_paged(api_key, sim, snap, mass_min_msun, N_max, dry_run=False):
    snap_data = api_get(f"{BASE_URL}/{sim}/snapshots/{snap}/", api_key)
    subs_url  = snap_data["subhalos"]
    log_min   = float(np.log10(mass_min_msun))

    params = {
        "limit":           100,
        "mass_stars__gte": log_min,
        "order_by":        "-mass_stars",
    }

    all_pos  = []
    page_url = subs_url
    print(f"  Fetching subhalos (paginated) for {sim} snap={snap}...")

    while page_url and len(all_pos) < N_max:
        if dry_run:
            print(f"  [DRY RUN] Would fetch page: {page_url}")
            break
        data     = api_get(page_url, api_key, params=params)
        for s in data.get("results", []):
            sd = api_get(s["url"], api_key)
            all_pos.append([sd["pos_x"], sd["pos_y"], sd["pos_z"]])
        page_url = data.get("next")
        params   = {}
        print(f"  {len(all_pos):,} / {N_max}", end="\r")

    if dry_run:
        return None

    pos_mpc = np.array(all_pos[:N_max]) * POS_UNIT
    print(f"\n  Fetched {len(pos_mpc):,} subhalos")
    return pos_mpc


# ── Main ──────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Fetch TNG100/TNG300 → universal .bin",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    p.add_argument("--api-key",  required=True)
    p.add_argument("--sim",      default="TNG100-1",
                   choices=list(TNG_BOXES.keys()))
    p.add_argument("--snap",     type=int,   default=99)
    p.add_argument("--mass-min", type=float, default=1e9)
    p.add_argument("--N-max",    type=int,   default=100000)
    p.add_argument("--method",   default="bulk", choices=["bulk","paged"])
    p.add_argument("--out",      default=None)
    p.add_argument("--inspect",  default=None)
    p.add_argument("--dry-run",  action="store_true")
    p.add_argument("--seed",     type=int, default=42)
    args = p.parse_args()

    if args.inspect:
        inspect_catalog(args.inspect)
        return

    np.random.seed(args.seed)

    if args.out is None:
        args.out = (f"../quantum_power_spectrum/data/"
                    f"{args.sim}_snap{args.snap}.bin")

    if not args.dry_run:
        print(f"Querying {args.sim} snapshot {args.snap}...")
        z_snap = get_snapshot_redshift(args.api_key, args.sim, args.snap)
        print(f"  z = {z_snap:.4f}")
    else:
        z_snap = 0.0
        print(f"[DRY RUN] {args.sim} snap={args.snap}")

    if args.method == "bulk":
        positions = fetch_bulk(
            args.api_key, args.sim, args.snap,
            args.mass_min, args.N_max, args.dry_run)
    else:
        positions = fetch_paged(
            args.api_key, args.sim, args.snap,
            args.mass_min, args.N_max, args.dry_run)

    if args.dry_run:
        print("\n[DRY RUN] Remove --dry-run to execute.")
        return

    weights = np.ones(len(positions))

    write_catalog(
        filename    = args.out,
        positions   = positions,
        weights     = weights,
        boxsize     = TNG_BOXES[args.sim],
        redshift    = z_snap,
        is_periodic = True,         # uniform periodic box
        has_weights = False,
        source_name = f"{args.sim}_snap{args.snap}",
        coord_type  = "cartesian_comoving",
    )
    inspect_catalog(args.out)


if __name__ == "__main__":
    main()
