"""
make_figure8.py -- ambiguity of a genuine simplex composition (K = 3).

Panel (a) controlled composition lattice, K = 3, L = 3 (10 states): the sharp
          interval for the exceedance probability P(Simpson >= 0.4) as a
          function of the number m of aggregated pixels.  The interval does not
          shrink with m; it widens.
Panel (b) the same lattice, the sharp interval for the mean Simpson diversity
          itself, with the independent and comonotonic couplings marked.
Panel (c) real ESA WorldCover, K = 3: block Simpson diversity against m, with
          the sharp band, the enumerated truth, and the independent prediction.

Reads results/compositional_controlled.csv and
      results/compositional_real_L3.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 120,
    "savefig.dpi": 320,
    "savefig.bbox": "tight",
    "legend.frameon": False,
    "legend.fontsize": 7.4,
})

C_BAND = "#8FB4D9"
C_LINE = "#2F5D8A"
C_TRUTH = "#B03A5B"
C_IND = "#C4762A"
C_COM = "#3E7C59"
C_REF = "#7A7873"


def fig8() -> None:
    ctrl = pd.read_csv(RES / "compositional_controlled.csv")
    real = pd.read_csv(RES / "compositional_real_L3.csv")
    ctrl = ctrl[ctrl["K"] == 3]
    real = real[real["levels"] == 3].sort_values("m")

    fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.5))

    # ---------------- (a) exceedance interval widens with m -------------
    ax = axs[0]
    ex = ctrl[ctrl["functional"] == "P(simpson>=0.4)"].sort_values("m")
    ms = ex["m"].to_numpy()
    lo = ex["mean_lower"].to_numpy()
    hi = ex["mean_upper"].to_numpy()
    for x, a, b in zip(ms, lo, hi):
        ax.plot([x, x], [a, b], "-", color=C_LINE, lw=6.0, alpha=0.30,
                solid_capstyle="round")
        ax.plot([x], [a], "o", ms=3.4, color=C_LINE)
        ax.plot([x], [b], "o", ms=3.4, color=C_LINE)
        ax.text(x, b + 0.035, f"{b - a:.2f}", ha="center", fontsize=6.8,
                color=C_REF)
    ax.set_xticks(ms)
    ax.set_xlabel("pixels aggregated, $m$")
    ax.set_ylabel(r"sharp interval for $\Pr(\mathrm{Simpson}\geq 0.4)$")
    ax.set_ylim(-0.05, 1.12)
    ax.set_title("(a) A nonlinear simplex query\nresolves less, not more")

    # ---------------- (b) mean Simpson interval -----------------------
    ax = axs[1]
    si = ctrl[ctrl["functional"] == "simpson"].sort_values("m")
    ms = si["m"].to_numpy()
    lo = si["mean_lower"].to_numpy()
    hi = si["mean_upper"].to_numpy()
    ax.fill_between(ms, lo, hi, color=C_BAND, alpha=0.55, lw=0,
                    label="sharp interval")
    ax.plot(ms, lo, "-", color=C_LINE, lw=1.3)
    ax.plot(ms, hi, "-", color=C_LINE, lw=1.3)
    ax.plot(si["m"], si["independent"], "-s", color=C_IND, ms=3.6, lw=1.2,
            label="independent coupling")
    ax.plot(si["m"], si["comonotonic"], "--", color=C_COM, lw=1.2,
            label="comonotonic coupling")
    ax.set_xticks(ms)
    ax.set_xlabel("pixels aggregated, $m$")
    ax.set_ylabel("mean Simpson diversity of a block")
    ax.set_ylim(0.24, 0.42)
    ax.legend(loc="upper right", ncol=1, frameon=True, facecolor="white",
              framealpha=0.92, edgecolor="none")
    ax.set_title("(b) The interval is set by information,\nnot by block size")

    # ---------------- (c) real WorldCover K = 3 ------------------------
    ax = axs[2]
    ms = real["m"].to_numpy()
    lo = real["lower"].to_numpy()
    hi = real["upper"].to_numpy()
    ax.fill_between(ms, lo, hi, color=C_BAND, alpha=0.55, lw=0,
                    label="sharp interval")
    ax.plot(ms, lo, "-", color=C_LINE, lw=1.3)
    ax.plot(ms, hi, "-", color=C_LINE, lw=1.3)
    ax.plot(ms, real["truth"], "o-", color=C_TRUTH, ms=4.6, lw=1.8, zorder=5,
            label="enumerated truth")
    ax.plot(ms, real["independent"], "-^", color=C_IND, ms=4.0, lw=1.4,
            label="independent coupling")
    ax.plot(ms, real["comonotonic"], "--", color=C_COM, lw=1.2,
            label="comonotonic coupling")
    ax.set_xticks(ms)
    ax.set_xlabel("pixels aggregated ($30$ m pixels), $m$")
    ax.set_ylabel("Simpson diversity of a block")
    ax.set_ylim(-0.03, 0.78)
    ax.set_title("(c) ESA WorldCover, $K=3$:\ntruth inside, independence far out")
    ax.legend(loc="upper left", ncol=2, columnspacing=0.8, handlelength=1.4,
              handletextpad=0.5, frameon=True, facecolor="white",
              framealpha=0.92, edgecolor="none")
    for x, t, i in zip(ms, real["truth"], real["independent"]):
        ax.annotate(f"{i / t:.1f}$\\times$", (x, i), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=6.6, color=C_IND)

    fig.tight_layout()
    fig.savefig(FIG / "fig8_compositional.png")
    plt.close(fig)
    print("wrote fig8_compositional.png")


if __name__ == "__main__":
    fig8()
