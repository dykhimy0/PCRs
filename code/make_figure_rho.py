"""make_figure_rho.py -- Figure 11: a single stored dependence scalar does not
transfer across block size.

Section 6.7 says that one exchangeable correlation rho turns the identified
interval into a point.  That is true, but it is only half the contract, and the
half that was missing in v9.  Exchangeability means the *same* rho reproduces
the block variance at every block size, so a producer who stores one number is
entitled to have it used at every block size.  On the eleven reference windows
it is not: the rho that reproduces the enumerated truth at m = 4 is not the rho
that reproduces it at m = 25, and the gap grows with m.

  (a) rho_exact(m) against m, one line per window.  Exchangeability would make
      every line flat; instead every line falls, and falls monotonically.
  (b) the block variance implied by reusing the m = 4 scalar at m = 25, as a
      percentage error against the enumerated truth -- and, for comparison,
      the same quantity built from the natural lag-1 autocorrelation.
  (c) what that error costs a decision: the sample size implied by the lag-1
      scalar as a multiple of the true requirement.

Both candidate scalars lie inside the PSD-feasible range in every window, so
the defect is bias, not feasibility; the figure is drawn to make that visible
by plotting the feasible lower bound rather than by asserting it.

Reads results_v2/rho_scalar.csv.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from figstyle import apply_style, save  # noqa: E402
from windows_v2 import BIOME_ORDER  # noqa: E402

apply_style()

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results_v2"

BIOME_COLOR = {
    "temperate": "#2F5D8A", "cool-temperate": "#3E7C8A",
    "mediterranean": "#B4762A", "boreal": "#5B4B8A",
    "alpine": "#7A6A4F", "subtropical": "#2E7D5B", "tropical": "#A8474B",
}
BIOME_SHORT = {
    "temperate": "temperature", "cool-temperate": "cool-temp",
    "mediterranean": "mediterr.", "boreal": "boreal", "alpine": "alpine",
    "subtropical": "subtrop.", "tropical": "tropical",
}

SHORT = {
    "N48E006 (26k,26k)": "Black Forest",
    "N48E006 (28k,20k)": "Kraichgau",
    "N48E009 (26k,26k)": "Swabian Jura",
    "N39W009 (18k,18k)": "Alentejo",
    "N36W006 (15k,15k)": "Andalusia",
    "N66E024 (26k,20k)": "Lapland",
    "N63E027 (14k,6k)": "Kainuu",
    "N21E099 (24k,12k)": "Xishuangbanna",
    "N24E117 (14k,6k)": "Wuyi Shan",
    "N27E093 (10k,16k)": "Bomi",
    "N45E126 (2k,25k)": "Xiaoxing'anling",
}

M_STORE, M_USE = 4, 25


def order_of(df: pd.DataFrame) -> list[str]:
    wins = df[["window", "biome"]].drop_duplicates()
    wins["k"] = wins["biome"].map(BIOME_ORDER.index)
    wins = wins.sort_values(["k", "window"])
    return list(wins["window"])


def main() -> None:
    df = pd.read_csv(RES / "rho_scalar.csv")
    order = order_of(df)
    m_star = df["m"].max()

    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.95),
                            gridspec_kw={"width_ratios": [1.0, 1.30, 1.0]})

    # ---------------- (a) rho_exact(m) is not constant --------------------
    ax = axs[0]
    for w in order:
        g = df[df["window"] == w].sort_values("m")
        col = BIOME_COLOR[g["biome"].iloc[0]]
        ax.plot(g["m"], g["rho_exact"], "-o", color=col, lw=1.1, ms=3.2,
                mec="white", mew=0.6, zorder=3)
    spread = df.groupby("window")["rho_exact_spread"].first()
    ax.annotate("exchangeability would\nmake every line flat",
                xy=(11.5, 0.80), xytext=(4.6, 0.955), fontsize=8.1,
                color="#5F5E5A", ha="left", va="top", zorder=4,
                arrowprops=dict(arrowstyle="->", color="#9A9A93", lw=0.9))
    ax.set_xticks([4, 9, 16, 25])
    ax.set_xticklabels(["4", "9", "16", "25"])
    ax.set_xlabel("pixels per block, $m$")
    ax.set_ylabel(r"scalar that reproduces the truth, $\rho_{\rm exact}$")
    ax.set_ylim(0.52, 0.98)
    ax.set_title("(a) the reproducing scalar falls with $m$", fontsize=9.6)
    ax.text(0.03, 0.045,
            f"spread {spread.min():.2f}\u2013{spread.max():.2f} across $m$",
            transform=ax.transAxes, fontsize=7.9, color="#5F5E5A",
            va="bottom", ha="left")

    # ---------------- (b) reusing the m = 4 scalar at m = 25 --------------
    ax = axs[1]
    rows = []
    for w in order:
        g = df[(df["window"] == w)]
        b = g["biome"].iloc[0]
        rho4 = float(g[g["m"] == M_STORE]["rho_exact"].iloc[0])
        sig2 = float(g["sigma2"].iloc[0])
        truth = float(g[g["m"] == M_USE]["truth_block_var"].iloc[0])
        v_app = sig2 * (1.0 + (M_USE - 1) * rho4) / M_USE
        rows.append(dict(window=w, biome=b,
                         err_store=100 * (v_app - truth) / truth,
                         err_lag1=float(g[g["m"] == M_USE]["err_lag1_pct"].iloc[0])))
    tb = pd.DataFrame(rows)
    x = np.arange(len(tb))
    cols = [BIOME_COLOR[b] for b in tb["biome"]]
    ax.bar(x - 0.19, tb["err_store"], width=0.36, color=cols, alpha=0.95,
           edgecolor="white", lw=0.5, label=f"scalar stored at $m={M_STORE}$")
    ax.bar(x + 0.19, tb["err_lag1"], width=0.36, color="none",
           edgecolor=cols, lw=1.1, hatch="///",
           label="lag-1 autocorrelation")
    ax.axhline(0, color="#6E6E68", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[w] for w in tb["window"]], rotation=45,
                       ha="right", fontsize=7.4)
    ax.set_ylabel("error in block variance at $m=25$ (%)")
    ax.set_title("(b) the stored scalar overstates the variance", fontsize=9.6)
    ax.legend(loc="upper left", fontsize=7.6, frameon=False)
    ax.set_ylim(0, max(tb["err_lag1"].max(), tb["err_store"].max()) * 1.32)
    ax.text(0.985, 0.955,
            f"+{tb['err_store'].min():.1f}% to +{tb['err_store'].max():.1f}%",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.9,
            color="#444441")

    # ---------------- (c) the sampling cost of the natural scalar ---------
    ax = axs[2]
    cr = np.array([float(df[(df["window"] == w) & (df["m"] == M_USE)]
                         ["cost_ratio"].iloc[0]) for w in order])
    ax.barh(np.arange(len(order)), cr, height=0.62,
            color=[BIOME_COLOR[df[df["window"] == w]["biome"].iloc[0]]
                   for w in order], alpha=0.95, edgecolor="white", lw=0.5)
    ax.axvline(1.0, color="#6E6E68", lw=0.9, ls="--")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([SHORT[w] for w in order], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("samples, lag-1 scalar $\\div$ truth")
    ax.set_title("(c) a plausible scalar costs extra samples", fontsize=9.6)
    ax.set_xlim(1.0, max(cr) * 1.06)
    ax.text(0.97, 0.05, f"{cr.min():.2f}\u2013{cr.max():.2f}\u00d7",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.9,
            color="#444441")

    fig.subplots_adjust(left=0.075, right=0.99, top=0.885, bottom=0.235,
                        wspace=0.42)
    save(fig, "fig_rho_transfer.png")
    print(f"  (a) rho_exact spread per window: "
          f"{spread.min():.4f} to {spread.max():.4f}")
    print(f"  (b) reuse m={M_STORE} at m={M_USE}: "
          f"+{tb['err_store'].min():.2f}% to +{tb['err_store'].max():.2f}% "
          f"(median +{tb['err_store'].median():.2f}%)")
    print(f"      lag-1 at m={M_USE}: +{tb['err_lag1'].min():.2f}% to "
          f"+{tb['err_lag1'].max():.2f}%")
    print(f"  (c) sample cost ratio: {cr.min():.3f} to {cr.max():.3f}")


if __name__ == "__main__":
    main()
