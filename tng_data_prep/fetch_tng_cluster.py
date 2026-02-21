#!/usr/bin/env python3
"""
fetch_tng_cluster.py
────────────────────
Download IllustrisTNG-Cluster halo catalog via the REST API and
write a universal .bin for the quantum_power_spectrum C++ project.

The correct approach for TNG-Cluster:
  - FoF halo positions   → GET /subhalos/?primary_flag=1  (central BCG = 1 per cluster)
  - Subhalo positions    → GET /subhalos/?primary_flag=0  (satellite subhalos)
  - FoF group fields     → GET /halos/N/info.json         (per-halo, gives GroupPos etc.)

NOTE: /halos/ and /groupcat/ bulk endpoints return HTML for TNG-Cluster.
      The only supported approach is the paginated /subhalos/ REST endpoint.

MODES:
  clusters  — one central subhalo per FoF halo (primary_flag=1)
  subhalos  — all satellite subhalos (primary_flag=0)
  combined  — clusters + subhalos with mass weights
"""

import argparse, os, sys, time
import numpy as np
import requests

from catalog_io import write_catalog, inspect_catalog

BASE_URL  = "https://www.tng-project.org/api"
SIM_NAME  = "TNG-Cluster"
TNG_H     = 0.6774
MASS_UNIT = 1e10 / TNG_H    # code mass → Msun
POS_UNIT  = 1.0  / 1000.0   # ckpc/h   → cMpc/h
TNG_BOX   = 680.0            # virtual container [cMpc/h]


# ── Core HTTP helper — JSON only ──────────────────────────────────────
def api_get(url, api_key, params=None, retries=4):
    headers = {"api-key": api_key}
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers,
                             params=params, timeout=120)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 20))
                print(f"  Rate limited — waiting {wait}s...")
                time.sleep(wait)
                continue
            r.raise_for_status()
            ct = r.headers.get("content-type", "")
            if "application/json" not in ct:
                raise RuntimeError(
                    f"Expected JSON, got: {ct}\n"
                    f"URL: {url}\n"
                    f"Body preview: {r.text[:400]}")
            return r.json()
        except requests.exceptions.Timeout:
            print(f"  Timeout attempt {attempt+1}/{retries}, retrying...")
            time.sleep(5 * (attempt + 1))
        except requests.exceptions.ConnectionError as e:
            print(f"  Connection error attempt {attempt+1}/{retries}: {e}")
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"API failed after {retries} attempts: {url}")


# ── Snapshot info ─────────────────────────────────────────────────────
def get_snapshot_info(api_key, snap):
    data = api_get(f"{BASE_URL}/{SIM_NAME}/snapshots/{snap}/", api_key)
    return {
        "redshift":  data["redshift"],
        "num_fof":   data.get("num_groups_fof",    0),
        "num_sub":   data.get("num_groups_subfind", 0),
        "subs_url":  data["subhalos"],   # ← comes directly from snapshot JSON
    }


# ── Paginated subhalo fetcher (core of both modes) ────────────────────
def _fetch_subhalos_paginated(api_key, subs_url, params, N_max,
                               label, dry_run=False):
    """
    Generic paginated fetcher for /subhalos/ endpoint.
    Returns list of (pos_x, pos_y, pos_z, mass) tuples in code units.
    """
    if dry_run:
        print(f"  [DRY RUN] GET {subs_url}")
        print(f"  [DRY RUN] params = {params}")
        return []

    all_records = []
    url         = subs_url
    page        = 0
    page_params = params.copy()

    while url and len(all_records) < N_max:
        data    = api_get(url, api_key, params=page_params)
        results = data.get("results", [])

        if page == 0:
            total = data.get("count", "?")
            print(f"  Total {label} matching filter: {total}")
            if total == 0:
                print(f"  WARNING: 0 results — check filter parameters")
                break

        for s in results:
            # Each result in the list has: id, mass_log_msun, url
            # We need the full subhalo record for pos_x/y/z and mass
            sub = api_get(s["url"], api_key)

            pos_x  = sub.get("pos_x")
            pos_y  = sub.get("pos_y")
            pos_z  = sub.get("pos_z")
            mass   = sub.get("mass")         # total mass in code units

            # Skip if any field is missing
            if None in (pos_x, pos_y, pos_z, mass):
                continue

            all_records.append((
                float(pos_x),
                float(pos_y),
                float(pos_z),
                float(mass)
            ))

            if len(all_records) >= N_max:
                break

        page      += 1
        url        = data.get("next")
        page_params = {}   # 'next' URL already contains all query params
        print(f"  Page {page}: {len(all_records):,} / {N_max} {label}",
              end="\r", flush=True)

    print(f"\n  Collected: {len(all_records):,} {label}")
    return all_records


# ── Mode 1: cluster positions (primary_flag=1 → one per FoF halo) ─────
def fetch_cluster_positions(api_key, snap_info, mass_min_msun,
                             N_max, dry_run=False):
    """
    primary_flag=1 → central (most massive) subhalo of each FoF group.
    One per cluster — this is the cluster position catalog.
    Filter by total mass >= mass_min.
    """
    mass_min_code = mass_min_msun / MASS_UNIT

    params = {
        "primary_flag": 1,
        "mass__gte":    mass_min_code,
        "order_by":     "-mass",
        "limit":        100,
    }

    print(f"  Fetching cluster centrals (primary_flag=1,"
          f" M > {mass_min_msun:.1e} Msun)...")

    records = _fetch_subhalos_paginated(
        api_key, snap_info["subs_url"], params, N_max,
        label="clusters", dry_run=dry_run)

    if dry_run or not records:
        return None, None

    pos_mpc = np.array([[r[0], r[1], r[2]] for r in records],
                       dtype=np.float64) * POS_UNIT
    masses  = np.array([r[3] * MASS_UNIT for r in records],
                       dtype=np.float64)
    return pos_mpc, masses


# ── Mode 2: satellite subhalo positions (primary_flag=0) ──────────────
def fetch_subhalo_positions(api_key, snap_info, sub_mass_min_msun,
                             N_max, dry_run=False):
    """
    primary_flag=0 → satellite subhalos (not the central).
    Filter by stellar mass >= sub_mass_min.
    """
    # API uses log10(mass_stars) in solar units for this filter
    sub_mass_min_log = np.log10(sub_mass_min_msun)

    params = {
        "primary_flag":    0,
        "mass_log_msun__gte": sub_mass_min_log,
        "order_by":        "-mass",
        "limit":           100,
    }

    print(f"  Fetching satellite subhalos (primary_flag=0,"
          f" M* > {sub_mass_min_msun:.1e} Msun)...")

    records = _fetch_subhalos_paginated(
        api_key, snap_info["subs_url"], params, N_max,
        label="subhalos", dry_run=dry_run)

    if dry_run or not records:
        return None, None

    pos_mpc = np.array([[r[0], r[1], r[2]] for r in records],
                       dtype=np.float64) * POS_UNIT
    masses  = np.array([r[3] * MASS_UNIT for r in records],
                       dtype=np.float64)
    return pos_mpc, masses


# ── Mode 3: combined ──────────────────────────────────────────────────
def fetch_combined(api_key, snap_info, mass_min_msun, sub_mass_min_msun,
                   N_clusters, N_subs_max, dry_run=False):
    pos_cl, m_cl = fetch_cluster_positions(
        api_key, snap_info, mass_min_msun, N_clusters, dry_run)
    pos_sh, m_sh = fetch_subhalo_positions(
        api_key, snap_info, sub_mass_min_msun, N_subs_max, dry_run)

    if dry_run:
        return None, None

    positions = np.vstack([pos_cl, pos_sh])
    weights   = np.concatenate([
        np.log10(m_cl + 1.0),
        np.ones(len(pos_sh))
    ])
    print(f"\n  Combined: {len(pos_cl)} clusters"
          f" + {len(pos_sh)} subhalos = {len(positions)} total")
    return positions, weights


# ── Stats printer ─────────────────────────────────────────────────────
def print_stats(positions, masses=None, label="Catalog"):
    print(f"\n── {label} Statistics {'─'*28}")
    print(f"  N objects  : {len(positions):,}")
    for i, ax in enumerate(["x", "y", "z"]):
        lo, hi = positions[:, i].min(), positions[:, i].max()
        print(f"  {ax} range   : [{lo:.2f}, {hi:.2f}] cMpc/h")
    if masses is not None:
        print(f"  mass range : [{masses.min():.2e}, {masses.max():.2e}] Msun")
        print(f"  median M   : {np.median(masses):.2e} Msun")
    print(f"{'─'*50}\n")


# ── Main ──────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Fetch TNG-Cluster → universal .bin")

    p.add_argument("--api-key",      required=True)
    p.add_argument("--snap",         type=int,   default=99)
    p.add_argument("--mode",         default="clusters",
                   choices=["clusters", "subhalos", "combined"])
    p.add_argument("--mass-min",     type=float, default=1e14)
    p.add_argument("--sub-mass-min", type=float, default=1e9)
    p.add_argument("--N-clusters",   type=int,   default=352)
    p.add_argument("--N-subs-max",   type=int,   default=500)
    p.add_argument("--out",          default=None)
    p.add_argument("--inspect",      default=None)
    p.add_argument("--dry-run",      action="store_true")
    p.add_argument("--seed",         type=int,   default=42)
    args = p.parse_args()

    print(f"fetch_tng_cluster.py starting")
    print(f"  api-key : {args.api_key[:6]}... ({len(args.api_key)} chars)")
    print(f"  snap    : {args.snap}")
    print(f"  mode    : {args.mode}")
    print(f"  out     : {args.out}")
    sys.stdout.flush()

    if args.inspect:
        inspect_catalog(args.inspect)
        return

    np.random.seed(args.seed)

    if args.out is None:
        args.out = (f"../quantum_power_spectrum/data/"
                    f"TNG-Cluster_snap{args.snap}_{args.mode}.bin")

    # ── Snapshot metadata ─────────────────────────────────────────────
    print(f"\nQuerying {SIM_NAME} snapshot {args.snap}...")
    snap_info = get_snapshot_info(args.api_key, args.snap)
    z_snap    = snap_info["redshift"]
    print(f"  z = {z_snap:.4f}  |  FoF groups: {snap_info['num_fof']:,}"
          f"  |  Subfind: {snap_info['num_sub']:,}")
    print(f"  subhalos URL: {snap_info['subs_url']}")

    if args.dry_run:
        print(f"\n[DRY RUN] mode='{args.mode}' → {args.out}")
        if args.mode in ("clusters", "combined"):
            fetch_cluster_positions(
                args.api_key, snap_info, args.mass_min,
                args.N_clusters, dry_run=True)
        if args.mode in ("subhalos", "combined"):
            fetch_subhalo_positions(
                args.api_key, snap_info, args.sub_mass_min,
                args.N_subs_max, dry_run=True)
        print("\n[DRY RUN] Remove --dry-run to execute.")
        return

    # ── Dispatch ──────────────────────────────────────────────────────
    if args.mode == "clusters":
        positions, masses = fetch_cluster_positions(
            args.api_key, snap_info, args.mass_min, args.N_clusters)
        weights     = np.ones(len(positions))
        has_weights = False
        src_name    = f"TNG-Cluster_snap{args.snap}_clusters"

    elif args.mode == "subhalos":
        positions, masses = fetch_subhalo_positions(
            args.api_key, snap_info, args.sub_mass_min, args.N_subs_max)
        weights     = np.ones(len(positions))
        has_weights = False
        src_name    = f"TNG-Cluster_snap{args.snap}_subhalos"

    elif args.mode == "combined":
        positions, weights = fetch_combined(
            args.api_key, snap_info,
            args.mass_min, args.sub_mass_min,
            args.N_clusters, args.N_subs_max)
        masses      = None
        has_weights = True
        src_name    = f"TNG-Cluster_snap{args.snap}_combined"

    print_stats(
        positions,
        masses if args.mode != "combined" else None,
        label=args.mode)

    write_catalog(
        filename    = args.out,
        positions   = positions,
        weights     = weights,
        boxsize     = TNG_BOX,
        redshift    = z_snap,
        is_periodic = False,
        has_weights = has_weights,
        source_name = src_name,
        coord_type  = "cartesian_comoving",
    )
    inspect_catalog(args.out)


if __name__ == "__main__":
    main()
