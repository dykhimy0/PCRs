"""
make_figure_estimator_v12.py -- v12 split of the training-block figure.

The v11 figure combined two panels:

  (a) 90th-percentile relative error of the plug-in block variance vs K
  (b) interval contraction, band width / IQR

Section 6.4 in v12 keeps only the main-text conclusion (a few dozen blocks
beat the band) and moves the estimator design, the bias and IQR diagnostics
and the Xiaoxing'anling exception to Supplementary S1.  So the figure is split:

  fig12a_training_blocks.png -- panel (a) alone, main text, given the full
                                width so the eleven curves stay separable
  figS1_estimator_full.png   -- the original two-panel figure, Supplementary

Both read the same CSV.  No data change, no number change.
"""

from __future__ import annotations

import sys as _sys
from itertools import cycle
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "code"))
_sys.path.insert(0, str(_Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figstyle import apply_style, save

apply_style()

ROOT = _Path(__file__).resolve().parent.parent
RES = ROOT / "results_v2"

BIOME_COL = {
    "temperate": "#2F5D8A", "cool-temperate": "#3E7C8A",
    "mediterranean": "#B4762A", "boreal": "#5B4B8A",
    "alpine": "#7A6A4F", "subtropical": "#2E7D5B", "tropical": "#A8474B",
}


def window_styles(labels) -> dict:
    dashes = cycle(["-", "--", "-.", ":"])
    return {lab: next(dashes) for lab in labels}


def load():
    df = pd.read_csv(RES / "estimator_shrinkage.csv")
    unknown = sorted(set(df["biome"]) - set(BIOME_COL))
    if unknown:
        raise KeyError(f"no colour defined for biome(s): {unknown}")
    labels = sorted(df["window"].unique())
    return df, labels, window_styles(labels)


def panel_a(ax, df, labels, LS):
    """90th-percentile relative error of the plug-in block variance vs K."""
    for label in labels:
        g = df[(df["window"] == label) & (df["m"] == 16)]
        ax.plot(g["K"], g["err_p90_pct"], color=BIOME_COL[g["biome"].iloc[0]],
                linestyle=LS[label], marker="o", markersize=3.4, linewidth=1.4)
    ax.axhline(100, color="#888780", linestyle=":", linewidth=1.0)
    ax.axhline(86, color="#A8474B", linestyle="--", linewidth=1.0)
    # one annotation for both reference lines: at a 400 % axis the two lines sit
    # 14 units apart, too close for two separate labels at this print size
    ax.text(0.985, 0.325,
            "level-A band width (100 %)\nindependence error, mean (86 %)",
            color="#555555", fontsize=8.4, linespacing=1.45,
            transform=ax.transAxes, ha="right", va="bottom")
    ax.set_xscale("log")
    # 400 %, not 200 %: the Xiaoxing'anling window reaches 384 % at K = 8, and
    # that outlier is the point of the panel, so it must not be clipped
    ax.set_ylim(0, 400)
    ax.set_yticks([0, 100, 200, 300, 400])
    ax.set_xlabel("training blocks  $K$")
    ax.set_ylabel("90th-pct relative error (%)")
    ax.grid(which="both", axis="x", color="#DDDDDD", linewidth=0.5)


def panel_b(ax, df, labels, LS):
    """Interval contraction: level-A band width divided by the IQR."""
    for label in labels:
        g = df[(df["window"] == label) & (df["m"] == 16)]
        ax.plot(g["K"], g["shrinkage_iqr"], color=BIOME_COL[g["biome"].iloc[0]],
                linestyle=LS[label], marker="o", markersize=3.0, linewidth=1.3)
    ax.axhline(1, color="#888780", linestyle=":", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("training blocks  $K$")
    ax.set_ylabel("shrinkage (band width / IQR)")
    ax.grid(which="both", axis="x", color="#DDDDDD", linewidth=0.5)


def biome_legend(fig, y):
    handles = [Line2D([], [], color=c, linewidth=2.2, label=b)
               for b, c in BIOME_COL.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8.2,
               frameon=False, bbox_to_anchor=(0.5, y), handlelength=2.0,
               columnspacing=1.5, labelspacing=0.32)


def main():
    df, labels, LS = load()

    # --- main text: panel (a) alone, full width ---------------------------
    fig, ax = plt.subplots(figsize=(7.5, 3.35))
    panel_a(ax, df, labels, LS)
    biome_legend(fig, 0.005)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.975, bottom=0.245)
    save(fig, "fig12a_training_blocks.png")
    plt.close(fig)

    # --- supplementary: the original two panels ---------------------------
    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.30))
    panel_a(axs[0], df, labels, LS)
    axs[0].set_title("(a) a few blocks beat the band and independence",
                     fontsize=10.2)
    panel_b(axs[1], df, labels, LS)
    axs[1].set_title("(b) interval contraction,  $m = 16$", fontsize=10.2)
    biome_legend(fig, 0.002)
    fig.subplots_adjust(left=0.078, right=0.988, top=0.905, bottom=0.245,
                        wspace=0.26)
    save(fig, "figS1_estimator_full.png")
    plt.close(fig)

    p90 = df[(df["m"] == 16) & (df["K"] == 256)]["err_p90_pct"]
    shr = df[(df["m"] == 16) & (df["K"] == 256)]["shrinkage_iqr"]
    print(f"  K=256: p90 error {p90.min():.1f}-{p90.max():.1f}%, "
          f"shrinkage {shr.min():.1f}-{shr.max():.1f}x")
    print("wrote fig12a_training_blocks.png and figS1_estimator_full.png")


if __name__ == "__main__":
    main()
