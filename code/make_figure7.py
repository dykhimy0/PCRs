"""
make_figure7.py -- asymptotic sharpness of the PSD relaxation.

Separate from make_figures.py so that the figure can be regenerated without
touching the earlier figures.  All series are read from
results/psd_sweep_m.csv and results/psd_sweep_p.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from exp_psd import gap_closed, gap_envelope, rel_gap_envelope

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 120,
    "savefig.dpi": 320,
    "savefig.bbox": "tight",
    "legend.frameon": False,
    "legend.fontsize": 7.6,
})

C_ENV = "#5F5E5A"
C_05 = "#2F5D8A"
C_03 = "#B03A5B"
C_IRR = "#3E7C59"
C_REF = "#B4B2A9"


def fig7():
    df = pd.read_csv(RES / "psd_sweep_m.csv")
    fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.4))

    # (a) gap vs m, log-log, with the 1/(4 m^2) envelope
    ax = axs[0]
    ms = np.arange(2, 65)
    ax.plot(ms, gap_envelope(ms), "-", color=C_ENV, lw=1.8,
            label="envelope $1/(4m^2)$")
    for p, col, lab in ((0.5, C_05, "$p=0.5$"), (0.3, C_03, "$p=0.3$"),
                        (1 / np.pi, C_IRR, r"$p=1/\pi$")):
        g = df[df["p"].round(12) == round(p, 12)].sort_values("m")
        ax.plot(g["m"], g["gap"].clip(lower=1e-9), "-o", ms=2.4, lw=1.1,
                color=col, label=lab)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.set_xlabel("pixels aggregated, $m$")
    ax.set_ylabel("joint-mixability gap")
    ax.set_title("(a) The gap is bounded by $1/(4m^2)$")
    ax.legend(loc="lower left")

    # (b) relative gap vs m with reference slopes
    ax = axs[1]
    ax.plot(ms, rel_gap_envelope(0.5, ms), "-", color=C_05, lw=1.7,
            label="$1/m^2$  ($p=0.5$)")
    ax.plot(ms, rel_gap_envelope(0.3, ms), "-", color=C_03, lw=1.4,
            label="$1/(4m^2p(1-p))$  ($p=0.3$)")
    for p, col in ((0.5, C_05), (0.3, C_03)):
        g = df[(df["p"].round(12) == round(p, 12))].sort_values("m")
        ax.plot(g["m"], g["rel_gap"].clip(lower=1e-9), "o", ms=3.0,
                color=col, alpha=0.85)
    ax.axhline(0.01, color=C_REF, lw=1.0, ls=":")
    ax.text(2.2, 0.0115, "1 %", fontsize=7, color=C_ENV)
    ax.axvline(10, color=C_REF, lw=0.8, ls=":")
    ax.text(10.4, 0.35, "$m=10$", fontsize=7, color=C_ENV, rotation=90,
            va="top")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.set_ylim(1e-5, 0.5)
    ax.set_xlabel("pixels aggregated, $m$")
    ax.set_ylabel("relative gap, gap / $\\sigma^2$")
    ax.set_title("(b) Relative error falls as $m^{-2}$")
    ax.legend(loc="lower left")

    # (c) gap vs p at fixed m -- the two regimes
    ax = axs[2]
    ps = np.linspace(0.005, 0.995, 1200)
    cols = {4: "#1F4E79", 8: "#2E8B7A", 16: "#B03A5B"}
    for m in (4, 8, 16):
        g = np.array([gap_closed(float(p), m) for p in ps])
        ax.plot(ps, g, "-", lw=1.7, color=cols[m], label=f"$m={m}$")
        ax.axhline(gap_envelope(m), color=cols[m], lw=0.8, ls=":")
    for m in (4, 8, 16):
        lim = m * np.minimum(ps, 1 - ps)
        span = ps[lim > 1]
        if span.size:
            ax.axvspan(span.min(), span.max(), color="#F0EEE8", zorder=0)
    ax.set_xlabel("class proportion $p$")
    ax.set_ylabel("joint-mixability gap")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.0006, 0.0175)
    ax.set_title("(c) Zero until $m\\min(p,1-p)>1$, then bounded")
    ax.text(0.5, 0.0152, "shaded: PSD corner binds", fontsize=7,
            color="#6B6B6B", ha="center")
    ax.legend(loc="upper right", ncol=3, columnspacing=0.8, handlelength=1.2)
    fig.tight_layout()
    fig.savefig(FIG / "fig7_psd_asymptotic.png")
    plt.close(fig)
    print("wrote fig7_psd_asymptotic.png")


if __name__ == "__main__":
    fig7()
