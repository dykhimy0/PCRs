"""
make_figure9.py -- the geostatistical baseline on all six real windows.

The three original windows all sit in the temperate belt; two added windows sit
in a Mediterranean biome and one in a boreal biome, so the figure can show what
does and does not change when the biome changes.

Panel (a) the binned empirical variogram of every window, normalised by its own
          field variance, with the 60-pixel truncation window shaded.  The
          curves do not coincide: the boreal field decorrelates fastest and the
          Mediterranean one slowest, so a 60-pixel fit window truncates them by
          very different amounts.
Panel (b) error against the enumerated block variance for every estimator, over
          all 6 windows x 4 block sizes (24 cases).  symlog vertical axis so
          that both sign and magnitude are visible.
Panel (c) the decision consequence: inflation of the required sample size
          relative to independence, one line per window, with the sharp upper
          bound as a dotted band.
Panel (d) the same errors averaged within each biome.  This is the panel that
          makes the biome dependence visible: the ordering of the bars is the
          same for every variogram-based estimator and the ordering follows the
          fitted range, not the sample size.

Reads results/baseline_variogram_curves.csv, baseline_variogram.csv,
      baseline_decision.csv.
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
    "legend.frameon": True,
    "legend.fontsize": 7.0,
})

# window -> (biome, colour, linestyle)
STYLE = {
    "N48E006 (26k,26k)": ("temperate",     "#2F5D8A", "-"),
    "N48E006 (28k,20k)": ("temperate",     "#4E8A62", "--"),
    "N48E009 (26k,26k)": ("temperate",     "#8FB8D6", "-."),
    "N39W009 (18k,18k)": ("mediterranean", "#C4762A", "-"),
    "N36W006 (15k,15k)": ("mediterranean", "#B03A5B", "--"),
    "N66E024 (26k,20k)": ("boreal",        "#5B4B8A", "-"),
}
BIOME_COLOR = {"temperate": "#2F5D8A", "mediterranean": "#C4762A",
               "boreal": "#5B4B8A"}
BIOME_ORDER = ["temperate", "mediterranean", "boreal"]
BIOME_LABEL = {"temperate": "temperate (original 3)",
               "mediterranean": "Mediterranean (added 2)",
               "boreal": "boreal (added 1)"}
C_TRUNC = "#B4762A"
C_REF = "#8A8A83"
RMAX_TRUNC = 60

ESTIMATORS = [("err_truncated_pct", "60-px\nexp", "#C4762A"),
              ("err_wide_pct", "300-px\nexp", "#4E8A62"),
              ("err_nested_pct", "nested\nexp", "#7B5EA7"),
              ("err_empirical_pct", "empirical\nvariogram", "#2F5D8A"),
              ("err_independence_pct", "independence", C_REF)]


def fig9() -> None:
    cv = pd.read_csv(RES / "baseline_variogram_curves.csv")
    bv = pd.read_csv(RES / "baseline_variogram.csv")
    dec = pd.read_csv(RES / "baseline_decision.csv")
    windows = list(STYLE.keys())
    n_case = len(bv)

    fig, axs = plt.subplots(2, 2, figsize=(9.8, 7.0))

    # ---------------- (a) variograms, normalised ------------------------
    ax = axs[0, 0]
    ax.axvspan(0, RMAX_TRUNC, color="#F1EFE9", zorder=0)
    ax.axhline(1.0, color=C_REF, lw=1.0, ls="--")
    ax.text(112, 1.055, "sill $=\\sigma^2$", fontsize=6.8, color=C_REF)
    for w in windows:
        _b, col, ls = STYLE[w]
        g = cv[cv["window"] == w].sort_values("lag")
        ax.plot(g["lag"], g["gamma"] / g["sigma2"], ls, lw=1.6, color=col,
                label=w)
        reach = g[g["gamma"] / g["sigma2"] >= 0.95]
        if len(reach):
            x = float(reach["lag"].iloc[0])
            ax.plot([x], [0.95], "v", ms=5.0, color=col, zorder=5)
    ax.axvline(RMAX_TRUNC, color=C_TRUNC, lw=1.3, ls=":")
    ax.text(RMAX_TRUNC + 5, 0.10, "60-px fit window", fontsize=6.8,
            color=C_TRUNC, rotation=90, va="bottom")
    ax.text(215, 0.55, "triangles: first lag\nwith $\\gamma\\geq0.95\\sigma^2$",
            fontsize=6.3, color=C_REF, ha="center")
    ax.set_xlabel("lag (30 m pixels)")
    ax.set_ylabel("empirical $\\gamma(h)\\,/\\,\\sigma^2$")
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1.14)
    ax.legend(loc="lower right", framealpha=0.95, edgecolor="none", ncol=2)
    ax.set_title("(a) The sill is not reached at 60 px,\n"
                 "least so in the boreal window")

    # ---------------- (b) error of each estimator -----------------------
    ax = axs[0, 1]
    rng = np.random.default_rng(0)
    for i, (col, tick, colr) in enumerate(ESTIMATORS):
        v = bv[col].to_numpy()
        x = i + rng.uniform(-0.19, 0.19, size=v.size)
        ax.plot(x, v, "o", ms=3.5, color=colr, alpha=0.8)
        ax.plot([i - 0.36, i + 0.36], [np.median(v)] * 2, "-", lw=2.6,
                color=colr, zorder=6, solid_capstyle="butt")
    ax.axhline(0, color="#B4B2A9", lw=1.0)
    ax.set_yscale("symlog", linthresh=5)
    ax.set_yticks([-95, -50, -20, -5, 0, 5, 10])
    ax.set_yticklabels(["-95", "-50", "-20", "-5", "0", "+5", "+10"])
    ax.set_xticks(range(len(ESTIMATORS)))
    ax.set_xticklabels([e[1] for e in ESTIMATORS], fontsize=6.2)
    ax.set_xlim(-0.6, len(ESTIMATORS) - 0.4)
    ax.set_ylim(-140, 34)
    ax.set_ylabel("error against truth (%)")
    ax.set_title(f"(b) All {n_case} cases "
                 "(6 windows $\\times$ 4 sizes);\nbar = median")

    # ---------------- (c) sample-size inflation -------------------------
    ax = axs[1, 0]
    for w in windows:
        _b, col, ls = STYLE[w]
        d = dec[dec["window"] == w].sort_values("m")
        infl = d["n_empirical"] / d["n_independence"]
        upper = d["n_upper"] / d["n_independence"]
        ax.fill_between(d["m"], infl, upper, color=col, alpha=0.13, lw=0)
        ax.plot(d["m"], infl, "o", ms=4.0, color=col, ls=ls, lw=1.5, label=w)
        ax.plot(d["m"], upper, ":", color=col, lw=1.0, alpha=0.8)
    ax.axhline(1.0, color=C_REF, lw=1.2, ls="-.")
    ax.text(4.4, 1.06, "independence", fontsize=6.8, color=C_REF)
    ax.set_yscale("log")
    ax.set_xticks([4, 9, 16, 25])
    ax.set_ylim(0.8, 70)
    ax.set_xlabel("pixels per block, $m$")
    ax.set_ylabel("samples needed, relative to independence")
    ax.legend(loc="lower right", framealpha=0.95, edgecolor="none", ncol=2)
    ax.text(4.35, 64.0, "solid: variogram estimate\ndotted: sharp upper bound",
            fontsize=6.3, color=C_REF, va="top", ha="left")
    lo = float((dec["n_empirical"] / dec["n_independence"]).min())
    hi = float((dec["n_empirical"] / dec["n_independence"]).max())
    ax.set_title(f"(c) A 1 % target costs {lo:.1f} to {hi:.1f} times\n"
                 "more samples than independence implies")

    # ---------------- (d) NEW: error averaged by biome ------------------
    ax = axs[1, 1]
    absd = bv.assign(**{c: bv[c].abs() for c, _t, _r in ESTIMATORS})
    n_est = 4                       # leave out the (off-scale) independence
    width = 0.19
    ys = np.arange(n_est)
    for j, biome in enumerate(BIOME_ORDER):
        sub = absd[absd["biome"] == biome]
        vals = [sub[ESTIMATORS[k][0]].mean() for k in range(n_est)]
        pos = ys + (j - 1) * width
        ax.barh(pos, vals, height=width * 0.92, color=BIOME_COLOR[biome],
                label=BIOME_LABEL[biome])
    ax.set_xscale("log")
    ax.set_xlim(1.5, 60)
    ax.set_xticks([2, 5, 10, 20, 40])
    ax.set_xticklabels(["2", "5", "10", "20", "40"])
    ax.set_yticks(ys)
    ax.set_yticklabels([ESTIMATORS[k][1].replace("\n", " ")
                        for k in range(n_est)], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("mean |error| against truth (%), log scale")
    ax.legend(loc="lower right", framealpha=0.95, edgecolor="none")
    lo_b = absd[absd["biome"] == "boreal"]["err_truncated_pct"].mean()
    hi_b = absd[absd["biome"] == "mediterranean"]["err_truncated_pct"].mean()
    ax.set_title("(d) The truncated window hurts least\n"
                 f"in the boreal biome ({lo_b:.0f} %) and most\n"
                 f"in the Mediterranean one ({hi_b:.0f} %)")

    fig.tight_layout()
    fig.savefig(FIG / "fig9_baseline_decision.png")
    plt.close(fig)
    print(f"wrote fig9_baseline_decision.png  ({n_case} cases, "
          f"{len(windows)} windows, "
          f"{bv['biome'].nunique()} biomes)")


if __name__ == "__main__":
    fig9()
