#!/usr/bin/env python3
"""
plot_results.py
───────────────
Visualize quantum_power_spectrum C++ pipeline outputs.
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import uniform_filter1d

plt.rcParams.update({
    "figure.dpi":        150,
    "font.family":       "serif",
    "font.size":         12,
    "axes.labelsize":    13,
    "axes.titlesize":    13,
    "legend.fontsize":   10,
    "xtick.direction":   "in",
    "ytick.direction":   "in",
    "xtick.top":         True,
    "ytick.right":       True,
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "grid.linestyle":    "--",
})

COLORS = {
    "theory_cl":    "#2166ac",
    "theory_qu":    "#d6604d",
    "data_cl":      "#1a1a1a",
    "data_qu":      "#4dac26",
    "quantum_delta":"#762a83",
}

# ── Loaders ────────────────────────────────────────────────────────────
def load(filename, usecols=None):
    try:
        data = np.loadtxt(filename, comments="#", usecols=usecols)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"  [WARN] Could not load {filename}: {e}")
        return None


def load_pk_theory(path="output/power_spectrum_theory.dat"):
    d = load(path, usecols=(0, 1))
    if d is None: return None, None
    return d[:, 0], d[:, 1]


def load_pk_quantum(path="output/power_spectrum_theory_quantum.dat"):
    """Returns k, delta_P, P_classical, P_quantum (col 0,1,2,3)"""
    d = load(path, usecols=(0, 1, 2, 3))
    if d is None: return None, None, None, None
    k      = d[:, 0]
    dP     = d[:, 1]   # delta_P
    Pk_cl  = d[:, 2]   # P_classical
    # BUG FIX: recompute P_quantum from P_cl + delta_P
    # (avoids any column-order confusion)
    Pk_qu  = Pk_cl + dP
    return k, dP, Pk_cl, Pk_qu


def load_pk_data(path="output/power_spectrum_data.dat"):
    d = load(path, usecols=(0, 1, 2, 3))
    if d is None: return None, None, None, None
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def load_pk_data_quantum(path="output/power_spectrum_data_quantum.dat"):
    d = load(path, usecols=(0, 1, 2, 3))
    if d is None: return None, None, None, None
    k     = d[:, 0]
    dP    = d[:, 1]
    Pk_cl = d[:, 2]
    Pk_qu = Pk_cl + dP   # same fix
    return k, dP, Pk_cl, Pk_qu


def load_xi_theory_cl(path="output/xi_theory_classical.dat"):
    d = load(path, usecols=(0, 1))
    if d is None: return None, None
    return d[:, 0], d[:, 1]


def load_xi_theory_qu(path="output/xi_theory_quantum.dat"):
    d = load(path, usecols=(0, 1))
    if d is None: return None, None
    r, xi = d[:, 0], d[:, 1]
    # BUG FIX: smooth the quantum xi(r) to suppress Fourier ringing
    # Use a 5-point moving average on the raw values before plotting
    xi_smooth = uniform_filter1d(xi, size=5)
    return r, xi_smooth


def load_xi_data_cl(path="output/xi_data_classical.dat"):
    d = load(path, usecols=(0, 1, 2, 3, 4, 5))
    if d is None: return None, None, None
    r, xi, err = d[:, 0], d[:, 1], d[:, 2]
    mask = d[:, 3] > 0
    return r[mask], xi[mask], err[mask]


def load_xi_data_qu(path="output/xi_data_quantum.dat"):
    d = load(path, usecols=(0, 1))
    if d is None: return None, None
    r, xi = d[:, 0], d[:, 1]
    xi_smooth = uniform_filter1d(xi, size=5)
    return r, xi_smooth


# ── Plot 1: P(k) ───────────────────────────────────────────────────────
def plot_pk(ax, out_dir):
    # Use P_classical from the quantum file (col 2) as the reference —
    # both classical and quantum are then on the same bias scaling (b²P)
    k_thq, dP, Pk_thq_cl, Pk_thq = load_pk_quantum()
    k_d,   Pk_d_raw, Pk_d, Nm    = load_pk_data()

    if k_thq is not None:
        ax.loglog(k_thq, Pk_thq_cl,
                  color=COLORS["theory_cl"], lw=2.0,
                  label=r"Theory — classical $P(k)$", zorder=3)
        ax.loglog(k_thq, Pk_thq,
                  color=COLORS["theory_qu"], lw=2.0, ls="--",
                  label=r"Theory — quantum $P(k)$", zorder=3)

    if k_d is not None and len(k_d) > 0:
        Nm_safe = np.where(Nm > 0, Nm, 1)
        err     = np.abs(Pk_d) / np.sqrt(Nm_safe)
        mask    = Pk_d > 0
        if np.any(mask):
            ax.errorbar(k_d[mask], Pk_d[mask],
                        yerr=err[mask],
                        fmt="o", ms=5, color=COLORS["data_cl"],
                        ecolor="gray", elinewidth=1.2, capsize=3,
                        label="TNG-Cluster data (shot subtracted)",
                        zorder=4)

    ax.set_xlabel(r"$k$ [$h$/Mpc]")
    ax.set_ylabel(r"$P(k)$ [(Mpc/$h$)$^3$]")
    ax.set_title(r"Matter Power Spectrum $P(k)$")
    ax.legend(loc="upper right")
    ax.set_xlim([5e-3, 1.0])


# ── Plot 2: ξ(r) ──────────────────────────────────────────────────────
def plot_xi(ax, out_dir):
    # Load classical xi from the same file as quantum xi
    # to ensure consistent bias
    r_cl,  xi_cl          = load_xi_theory_cl()
    r_qu,  xi_qu          = load_xi_theory_qu()
    r_d,   xi_d,  xi_err  = load_xi_data_cl()
    r_dq,  xi_dq          = load_xi_data_qu()

    if r_cl is not None:
        mask = xi_cl > 1e-6
        ax.semilogy(r_cl[mask], xi_cl[mask],
                    color=COLORS["theory_cl"], lw=2.0,
                    label=r"Theory — classical $\xi(r)$")

    if r_qu is not None:
        mask = xi_qu > 1e-6
        if np.any(mask):
            ax.semilogy(r_qu[mask], xi_qu[mask],
                        color=COLORS["theory_qu"], lw=2.0, ls="--",
                        label=r"Theory — quantum $\xi(r)$")

    if r_d is not None and len(r_d) > 0:
        ax.errorbar(r_d, xi_d,
                    yerr=xi_err,
                    fmt="o", ms=5, color=COLORS["data_cl"],
                    ecolor="gray", elinewidth=1.2, capsize=3,
                    label=r"TNG-Cluster — Landy-Szalay $\xi(r)$",
                    zorder=4)

    if r_dq is not None:
        mask = xi_dq > 1e-6
        if np.any(mask):
            ax.semilogy(r_dq[mask], xi_dq[mask],
                        color=COLORS["data_qu"], lw=1.5, ls="-.",
                        label=r"TNG-Cluster — quantum $\xi(r)$")

    ax.set_xlabel(r"$r$ [Mpc/$h$]")
    ax.set_ylabel(r"$\xi(r)$")
    ax.set_title(r"Two-Point Correlation Function $\xi(r)$")
    ax.set_xscale("log")
    ax.legend(loc="upper right")
    ax.set_xlim([1, 300])
    ax.set_ylim([1e-3, 1e4])


# ── Plot 3: ΔP(k)/P(k) ────────────────────────────────────────────────
def plot_delta_pk(ax, out_dir):
    k_thq, dP, Pk_cl, Pk_qu   = load_pk_quantum()
    k_dq,  dPq, Pk_dcl, Pk_dqu = load_pk_data_quantum()

    if k_thq is not None and np.any(Pk_cl > 0):
        safe = np.where(Pk_cl > 0, Pk_cl, 1.0)
        frac_th = dP / safe
        ax.semilogx(k_thq, frac_th * 100,
                    color=COLORS["theory_qu"], lw=2.0,
                    label=r"Theory: $\Delta P / P_{\rm cl}$")

    if k_dq is not None and np.any(Pk_dcl > 0):
        safe = np.where(Pk_dcl > 0, Pk_dcl, 1.0)
        frac_d = dPq / safe
        ax.semilogx(k_dq, frac_d * 100,
                    color=COLORS["data_qu"], lw=2.0, ls="--",
                    label=r"TNG-Cluster: $\Delta P / P_{\rm cl}$")

    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$k$ [$h$/Mpc]")
    ax.set_ylabel(r"$\Delta P(k)/P_{\rm classical}(k)$ [%]")
    ax.set_title(r"Fractional Quantum Correction to $P(k)$")
    ax.legend(loc="upper right")


# ── Plot 4: BAO — r²ξ(r) ─────────────────────────────────────────────
def plot_bao(ax, out_dir):
    r_cl,  xi_cl          = load_xi_theory_cl()
    r_qu,  xi_qu          = load_xi_theory_qu()
    r_d,   xi_d,  xi_err  = load_xi_data_cl()

    r_min, r_max = 50.0, 200.0

    if r_cl is not None:
        mask = (r_cl >= r_min) & (r_cl <= r_max)
        ax.plot(r_cl[mask], r_cl[mask]**2 * xi_cl[mask],
                color=COLORS["theory_cl"], lw=2.0,
                label="Theory — classical")

    if r_qu is not None:
        mask = (r_qu >= r_min) & (r_qu <= r_max)
        if np.any(mask):
            ax.plot(r_qu[mask], r_qu[mask]**2 * xi_qu[mask],
                    color=COLORS["theory_qu"], lw=2.0, ls="--",
                    label="Theory — quantum")

    if r_d is not None and len(r_d) > 0:
        mask = (r_d >= r_min) & (r_d <= r_max)
        if np.any(mask):
            ax.errorbar(r_d[mask], r_d[mask]**2 * xi_d[mask],
                        yerr=r_d[mask]**2 * xi_err[mask],
                        fmt="o", ms=5, color=COLORS["data_cl"],
                        ecolor="gray", elinewidth=1.2, capsize=3,
                        label="TNG-Cluster")

    ax.axvline(105.0, color="k", lw=0.8, ls=":",
               label=r"BAO peak $\sim$105 Mpc/$h$")
    ax.set_xlabel(r"$r$ [Mpc/$h$]")
    ax.set_ylabel(r"$r^2\xi(r)$ [Mpc$^2$/$h^2$]")
    ax.set_title(r"BAO Feature: $r^2\xi(r)$")
    ax.set_xlim([r_min, r_max])
    ax.legend(loc="upper left", fontsize=9)


# ── Summary panel ─────────────────────────────────────────────────────
def make_summary(out_dir):
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            hspace=0.35, wspace=0.30)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    plot_pk      (ax1, out_dir)
    plot_xi      (ax2, out_dir)
    plot_delta_pk(ax3, out_dir)
    plot_bao     (ax4, out_dir)

    fig.suptitle(
        "Quantum Field Theory Corrections to the Cosmic Power Spectrum\n"
        r"TNG-Cluster $z=0$, 352 clusters, $M_{200c} > 10^{14}\,M_\odot$",
        fontsize=13, y=1.01
    )
    path = os.path.join(out_dir, "summary_panel.png")
    fig.savefig(path, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Individual plots ──────────────────────────────────────────────────
def make_individual(out_dir):
    for name, fn in [
        ("power_spectrum",    plot_pk),
        ("correlation_function", plot_xi),
        ("quantum_correction",   plot_delta_pk),
        ("bao_feature",          plot_bao),
    ]:
        fig, ax = plt.subplots(figsize=(8, 6))
        fn(ax, out_dir)
        fig.tight_layout()
        path = os.path.join(out_dir, f"{name}.png")
        fig.savefig(path, bbox_inches="tight")
        print(f"  Saved: {path}")
        plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="output/plots")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    make_summary   (args.output_dir)
    make_individual(args.output_dir)
    print("\nDone.")
