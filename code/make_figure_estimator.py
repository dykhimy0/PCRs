"""
make_figure_estimator.py -- Figure: shrinking the identified interval from
training blocks (the concept proof of Section 6.8).

Reads results/estimator_shrinkage.csv (exp_estimator.py).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from figstyle import apply_style

apply_style()

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

BIOME_COL = {"temperate": "#1A6FB5", "mediterranean": "#D1823C", "boreal": "#3E7C59"}
WINDOW_LS = {"N48E006 (26k,26k)": "-", "N48E006 (28k,20k)": "--",
             "N48E009 (26k,26k)": "-.", "N39W009 (18k,18k)": "-",
             "N36W006 (15k,15k)": "--", "N66E024 (26k,20k)": "-"}


def main():
    df = pd.read_csv(RES / "estimator_shrinkage.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.9))

    # (a) 90th-percentile relative error of the plug-in block variance vs K
    ax = axs[0]
    for label, g in df.groupby("window"):
        g = g[g["m"] == 16]
        ax.plot(g["K"], g["err_p90_pct"], color=BIOME_COL[g["biome"].iloc[0]],
                linestyle=WINDOW_LS[label], marker="o", markersize=3.5, linewidth=1.4)
    ax.axhline(100, color="#888780", linestyle=":", linewidth=1.0)
    ax.text(4.2, 102, "level-A band width", color="#888780", fontsize=10, va="bottom")
    ax.axhline(86, color="#A8474B", linestyle="--", linewidth=1.0)
    ax.text(4.2, 88, "independence error (mean)", color="#A8474B", fontsize=10, va="bottom")
    ax.set_xscale("log")
    ax.set_ylim(0, 200)
    ax.set_xlabel("training blocks  K")
    ax.set_ylabel("90th-pct relative error (%)")
    ax.set_title("(a) a few blocks beat the band and independence", fontsize=12.5)
    ax.grid(which="both", axis="x", color="#DDDDDD", linewidth=0.5)

    # (b) shrinkage: level-A width divided by the IQR of the estimate
    ax = axs[1]
    for label, g in df.groupby("window"):
        g = g[g["m"] == 16]
        ax.plot(g["K"], g["shrinkage_iqr"], color=BIOME_COL[g["biome"].iloc[0]],
                linestyle=WINDOW_LS[label], marker="o", markersize=3.5, linewidth=1.4)
    ax.axhline(1, color="#888780", linestyle=":", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("training blocks  K")
    ax.set_ylabel("shrinkage (band width / IQR)")
    ax.set_title("(b) interval contraction,  m = 16", fontsize=12.5)
    ax.grid(which="both", axis="x", color="#DDDDDD", linewidth=0.5)

    # shared biome legend (draw once on panel a)
    handles = [plt.Line2D([], [], color=c, linewidth=2, label=b)
               for b, c in BIOME_COL.items()]
    axs[0].legend(handles=handles, loc="upper right", fontsize=10.5)

    fig.tight_layout()
    fig.savefig(FIG / "fig11_estimator.png")
    plt.close(fig)
    print("wrote fig11_estimator.png")


if __name__ == "__main__":
    main()
