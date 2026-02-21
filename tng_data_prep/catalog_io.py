"""
catalog_io.py — Universal catalog reader/writer.
All data preparation scripts import this module.

Binary format specification:
  HEADER (256 bytes):
    [0:8]    magic identifier  "QPSC_v1\0"
    [8:16]   int64   N_gal
    [16:24]  double  boxsize      (cMpc/h, 0.0 for surveys)
    [24:32]  double  redshift
    [32]     uint8   is_periodic  (1=simulation box, 0=survey)
    [33]     uint8   has_weights
    [34:98]  char64  source_name
    [98:162] char64  coord_type   ("cartesian_comoving" | "radec_redshift")
    [162:256] reserved (zeros)
  DATA (N_gal × 32 bytes):
    double x, double y, double z, double weight
"""
import struct
import os
import numpy as np

HEADER_SIZE   = 256
RECORD_SIZE   = 32        # 4 × 8 bytes
MAGIC         = b'QPSC_v1\x00'


def write_catalog(filename: str,
                  positions:    np.ndarray,   # (N, 3)
                  weights:      np.ndarray,   # (N,)
                  boxsize:      float,
                  redshift:     float,
                  is_periodic:  bool,
                  has_weights:  bool,
                  source_name:  str,
                  coord_type:   str = "cartesian_comoving"):
    """Write a universal binary catalog readable by the C++ project."""

    N = len(positions)
    assert positions.shape == (N, 3), "positions must be (N, 3)"
    assert weights.shape   == (N,),   "weights must be (N,)"

    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)

    with open(filename, "wb") as f:

        # ── Build header ─────────────────────────────────────────────
        header = bytearray(HEADER_SIZE)

        header[0:8]   = MAGIC
        struct.pack_into("q", header,  8, N)
        struct.pack_into("d", header, 16, boxsize)
        struct.pack_into("d", header, 24, redshift)
        struct.pack_into("B", header, 32, int(is_periodic))
        struct.pack_into("B", header, 33, int(has_weights))

        sname = source_name.encode("utf-8")[:64].ljust(64, b'\x00')
        header[34:98]   = sname

        ctype = coord_type.encode("utf-8")[:64].ljust(64, b'\x00')
        header[98:162]  = ctype
        # bytes 162–255 stay zero (reserved)

        f.write(header)

        # ── Write data as raw float64 array ──────────────────────────
        data = np.column_stack([positions, weights]).astype(np.float64)
        data.tofile(f)

    size_mb = os.path.getsize(filename) / 1e6
    print(f"\nWritten : {filename}  ({N:,} galaxies,  {size_mb:.2f} MB)")
    print(f"  source_name  : {source_name}")
    print(f"  coord_type   : {coord_type}")
    print(f"  is_periodic  : {is_periodic}")
    print(f"  boxsize      : {boxsize} cMpc/h")
    print(f"  redshift     : {redshift:.6f}")


def read_catalog_header(filename: str) -> dict:
    """Read only the header without loading all data."""
    with open(filename, "rb") as f:
        raw = f.read(HEADER_SIZE)

    magic    = raw[0:8]
    if magic != MAGIC:
        raise ValueError(f"Bad magic in {filename}: {magic}")

    N        = struct.unpack_from("q", raw,  8)[0]
    boxsize  = struct.unpack_from("d", raw, 16)[0]
    redshift = struct.unpack_from("d", raw, 24)[0]
    is_per   = struct.unpack_from("B", raw, 32)[0]
    has_w    = struct.unpack_from("B", raw, 33)[0]
    src_name = raw[34:98].rstrip(b'\x00').decode("utf-8")
    ctype    = raw[98:162].rstrip(b'\x00').decode("utf-8")

    return {
        "N_gal":       N,
        "boxsize":     boxsize,
        "redshift":    redshift,
        "is_periodic": bool(is_per),
        "has_weights": bool(has_w),
        "source_name": src_name,
        "coord_type":  ctype,
    }


def read_catalog_data(filename: str):
    """
    Read full catalog. Returns (positions, weights, header_dict).
    positions shape: (N, 3) float64
    weights   shape: (N,)   float64
    """
    header = read_catalog_header(filename)
    N      = header["N_gal"]

    with open(filename, "rb") as f:
        f.seek(HEADER_SIZE)
        raw  = np.frombuffer(f.read(N * RECORD_SIZE), dtype=np.float64)

    data      = raw.reshape(N, 4)
    positions = data[:, :3]
    weights   = data[:,  3]

    return positions, weights, header


def inspect_catalog(filename: str):
    """Print a human-readable summary of a .bin catalog."""
    if not os.path.exists(filename):
        print(f"[ERROR] File not found: {filename}")
        return

    header = read_catalog_header(filename)
    size_mb = os.path.getsize(filename) / 1e6

    print(f"\n{'═'*50}")
    print(f"  Catalog : {os.path.basename(filename)}")
    print(f"{'═'*50}")
    print(f"  File size    : {size_mb:.2f} MB")
    for k, v in header.items():
        print(f"  {k:<14}: {v}")

    # Quick position stats
    positions, weights, _ = read_catalog_data(filename)
    for i, axis in enumerate(["x", "y", "z"]):
        lo, hi = positions[:, i].min(), positions[:, i].max()
        print(f"  {axis} range      : [{lo:.2f}, {hi:.2f}] cMpc/h")
    print(f"  weight range : [{weights.min():.4f}, {weights.max():.4f}]")
    print(f"{'═'*50}\n")
