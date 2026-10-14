# -*- coding: utf-8 -*-
"""
make_figure15_16.py -- the geostatistical route on the metric support.

Figure 15 (four panels)
  (a) normalised empirical variograms for the eleven windows, in metres, with
      the per-window truncated fit window shaded.  The x axis is a ground
      distance; on the geographic grid an integer lag mixed two distances that
      differ by 1/cos(lat), so this curve was not previously a distance at all.
  (b) error against the enumerated reference value for five estimators over all
      44 cases (11 windows x 4 block sizes).  symlog so sign and size are both
      visible; the bar is the median.
  (c) the decision consequence: pixels needed for a 1 percentage-point target,
      relative to independence.
  (d) mean absolute error by biome, for the truncated fit and the empirical
      curve.  The ordering follows the resolved fraction of the sill, not the
      sample size.

Figure 16 (two panels)
  (a) mean absolute bias of the truncated fit against the resolved fraction of
      the sill;  (b) the resolved fraction against the fitted range in metres.

Terminology: the reference quantity on WorldCover is the **enumerated reference
value** (the block value of the reference field, obtained by exhaustive
enumeration), never "truth" or "ground truth".  In the synthetic experiments of
Section 6.2 the constructed, error-free value is genuinely "truth"; the two must
not share a name.

Style: Times New Roman, STIX mathtext, large near-square panels, house type
sizes.  Drawings use the Python/matplotlib backend only.

Reads results_metric/metric_variogram_curves.csv, metric_variogram.csv,
      metric_decision.csv, metric_fits.csv, metric_blocks.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
sys.stdout.reconfigure(encoding="utf-8")

import matplotlib as mpl                              # noqa: E402
import matplotlib.pyplot as plt                       # noqa: E402

from figstyle import apply_style                      # noqa: E402

sys.path.insert(0, str(Path("C:/Users/ASUS/.workbuddy/skills/nature-figure/scripts")))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

RES = ROOT / "results_metric"
OUT_DIRS = (ROOT / "figures_v10", ROOT / "figures", ROOT / "figures_v2")

# --- house type sizes, matched to Figure 13 -------------------------------
# Raised to the same scale as make_figure13.py so that the whole figure set
# carries one type size: the panels are wide enough for the larger type, and
# the figure must stay readable once scaled into the journal column.
FS_TICK = 20.0
FS_LABEL = 23.0
FS_TITLE = 22.0
FS_LEGEND = 19.0
FS_ANNOT = 19.0

BIOME_COLOUR = {
    "boreal": "#3d6b8c", "cool-temperate": "#5b8c6b",
    "temperate": "#7ba05b", "mediterranean": "#c9a227",
    "alpine": "#8c7b9c", "subtropical": "#c2703d", "tropical": "#b5453a",
}
BIOME_ORDER = ["boreal", "cool-temperate", "temperate", "mediterranean",
               "alpine", "subtropical", "tropical"]
C_REF = "#8A8A83"

ESTIMATORS = [
    ("err_truncated_pct", "trunc.\nexp.", "#C4762A"),
    ("err_wide_pct", "wide\nexp.", "#4E8A62"),
    ("err_nested_pct", "nested\nexp.", "#7B5EA7"),
    ("err_empirical_pct", "empirical\nvariog.", "#2F5D8A"),
    ("err_independence_pct", "indep.", C_REF),
]
SHORT = {"boreal": "boreal", "cool-temperate": "cool-temp",
         "temperate": "temperate", "mediterranean": "mediterr.",
         "alpine": "alpine", "subtropical": "subtrop.",
         "tropical": "tropical"}


def _style(ax) -> None:
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(FS_TICK)
    ax.xaxis.label.set_fontsize(FS_LABEL)
    ax.yaxis.label.set_fontsize(FS_LABEL)
    ax.title.set_fontsize(FS_TITLE)


# ==========================================================================
# Figure 15
# ==========================================================================
def figure15() -> None:
    cur = pd.read_csv(RES / "metric_variogram_curves.csv")
    var = pd.read_csv(RES / "metric_variogram.csv")
    dec = pd.read_csv(RES / "metric_decision.csv")
    fits = pd.read_csv(RES / "metric_fits.csv")
    wins = sorted(cur["window"].unique())
    n_case = len(var)

    fig, axs = plt.subplots(2, 2, figsize=(15.6, 13.6),
                            gridspec_kw={"wspace": 0.26, "hspace": 0.46})

    # ---- (a) normalised variograms, metres -----------------------------
    ax = axs[0, 0]
    for w in wins:
        g = cur[cur["window"] == w].sort_values("lag_m")
        biome = g["biome"].iloc[0]
        ax.plot(g["lag_m"], g["gamma"] / g["sigma2"], lw=1.9,
                color=BIOME_COLOUR[biome], alpha=0.85)
    ax.axhline(1.0, color=C_REF, lw=1.2, ls="--")
    ax.text(7000, 1.06, "sill $=\\sigma^{2}$", fontsize=FS_ANNOT - 1.5,
            color=C_REF, ha="right")
    # shading: the truncated window is per window, so shade the median one
    med_trunc = float(fits.drop_duplicates("window")["rmax_m"].median()) * 0.20
    ax.axvspan(0, med_trunc, color="#EFEBE0", zorder=0)
    ax.text(med_trunc * 1.06, 0.10, "truncated\nfit window",
            fontsize=FS_ANNOT - 2.0, color="#8A6A2A", va="bottom")
    ax.set_xlabel("separation (m)")
    ax.set_ylabel(r"empirical $\gamma(h)\,/\,\sigma^{2}$")
    ax.set_xlim(0, 9000)
    ax.set_ylim(0, 1.22)
    ax.set_title("(a) The variogram is a function of a\n"
                 "ground distance, not of a pixel count", loc="left")
    _style(ax)

    # ---- (b) estimator errors over 44 cases ----------------------------
    ax = axs[0, 1]
    rng = np.random.default_rng(0)
    for i, (col, tick, colr) in enumerate(ESTIMATORS):
        v = var[col].to_numpy()
        x = i + rng.uniform(-0.20, 0.20, size=v.size)
        ax.plot(x, v, "o", ms=5.0, color=colr, alpha=0.72,
                markeredgewidth=0)
        ax.plot([i - 0.38, i + 0.38], [np.median(v)] * 2, "-", lw=3.2,
                color=colr, zorder=6, solid_capstyle="butt")
    ax.axhline(0, color="#B4B2A9", lw=1.2)
    ax.set_yscale("symlog", linthresh=5)
    ax.set_yticks([-95, -50, -20, -5, 0, 10])
    ax.set_yticklabels(["-95", "-50", "-20", "-5", "0", "+10"])
    ax.set_xticks(range(len(ESTIMATORS)))
    ax.set_xticklabels([e[1] for e in ESTIMATORS], fontsize=FS_TICK - 3.0)
    ax.set_xlim(-0.6, len(ESTIMATORS) - 0.4)
    ax.set_ylim(-150, 45)
    # label the median of the truncated fit and of independence, the two that matter
    ax.annotate(f"{np.median(var['err_truncated_pct']):.0f} %", xy=(0, -23),
                fontsize=FS_ANNOT - 1.0, color="#C4762A", ha="center",
                va="bottom")
    ax.annotate(f"{np.median(var['err_empirical_pct']):.1f} %", xy=(3, -2.8),
                fontsize=FS_ANNOT - 1.0, color="#2F5D8A", ha="center",
                va="top")
    ax.set_ylabel("error against the enumerated\nreference value (%)")
    ax.set_title(f"(b) All {n_case} cases "
                 f"({len(wins)} windows $\\times$ 4 sizes);\nbar = median",
                 loc="left")
    _style(ax)

    # ---- (c) sample-size inflation -------------------------------------
    ax = axs[1, 0]
    for w in wins:
        d = dec[dec["window"] == w].sort_values("m")
        col = BIOME_COLOUR[d["biome"].iloc[0]]
        infl = d["n_empirical"] / d["n_independence"]
        upper = d["n_upper"] / d["n_independence"]
        ax.fill_between(d["m"], infl, upper, color=col, alpha=0.12, lw=0)
        ax.plot(d["m"], infl, "o-", ms=6.0, lw=1.8, color=col)
        ax.plot(d["m"], upper, ":", color=col, lw=1.3, alpha=0.9)
    ax.axhline(1.0, color=C_REF, lw=1.4, ls="-.")
    ax.text(4.3, 1.08, "independence", fontsize=FS_ANNOT - 1.5, color=C_REF)
    ax.set_yscale("log")
    ax.set_xticks([4, 9, 16, 25])
    ax.set_xlim(3.2, 26.5)
    ax.set_ylim(0.8, 40)
    ax.set_xlabel("pixels per block, $m$")
    ax.set_ylabel("pixels needed, relative\nto independence")
    lo = float((dec["n_empirical"] / dec["n_independence"]).min())
    hi = float((dec["n_empirical"] / dec["n_independence"]).max())
    ax.set_title(f"(c) A 1 pp target costs {lo:.1f} to {hi:.1f} times\n"
                 "more pixels than independence implies", loc="left")
    _style(ax)

    # ---- (d) mean absolute error by biome -------------------------------
    ax = axs[1, 1]
    absd = var.assign(**{c: var[c].abs() for c, _t, _r in ESTIMATORS})
    biomes = [b for b in BIOME_ORDER if b in set(absd["biome"])]
    # order by resolved fraction of the sill (worst-resolved first)
    res_frac = (fits.drop_duplicates("window")
                .groupby("biome")["gamma_at_trunc_over_sill"].mean())
    biomes = sorted(biomes, key=lambda b: res_frac[b])
    ys = np.arange(len(biomes))
    width = 0.36
    for j, (key, lab, colr) in enumerate(
            [("err_truncated_pct", "truncated\nexponential", "#C4762A"),
             ("err_empirical_pct", "empirical\nvariogram", "#2F5D8A")]):
        vals = [absd[absd["biome"] == b][key].mean() for b in biomes]
        ax.barh(ys + (j - 0.5) * width, vals, height=width * 0.92,
                color=colr, label=lab)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{SHORT[b]}  ({res_frac[b]:.2f})" for b in biomes],
                       fontsize=FS_TICK - 3.5)
    ax.invert_yaxis()
    ax.set_xlabel("mean |error| (%), log scale")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 90)
    ax.set_xticks([1, 3, 10, 30])
    ax.set_xticklabels(["1", "3", "10", "30"])
    ax.legend(loc="lower right", fontsize=FS_LEGEND)
    ax.set_title("(d) The bias is graded by the resolved\n"
                 "fraction in brackets, not by biome", loc="left")
    _style(ax)

    fig.subplots_adjust(left=0.098, right=0.985, top=0.925, bottom=0.088)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=ROOT / "figures_v10" / "fig15.panel-alignment.json",
        overlay_svg=ROOT / "figures_v10" / "fig15.panel-alignment.svg",
    )
    for d in OUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / "fig15_variogram_baseline.png")
        fig.savefig(d / "fig15_variogram_baseline.pdf")
        fig.savefig(d / "fig15_variogram_baseline.svg")
    plt.close(fig)
    print("wrote fig15_variogram_baseline.png/.pdf/.svg")


# ==========================================================================
# Figure 16
# ==========================================================================
def figure16() -> None:
    fits = pd.read_csv(RES / "metric_fits.csv").drop_duplicates("window")
    var = pd.read_csv(RES / "metric_variogram.csv")

    # mean absolute truncated-fit bias per window (over the four block sizes)
    bias = (var.assign(a=var["err_truncated_pct"].abs())
            .groupby("window", as_index=False)["a"].mean())
    m = fits.merge(bias, on="window")
    m["resolved"] = m["gamma_at_trunc_over_sill"]

    fig, axs = plt.subplots(1, 2, figsize=(16.4, 8.6),
                            gridspec_kw={"wspace": 0.28})

    # ---- (a) bias vs resolved fraction --------------------------------
    ax = axs[0]
    for _, r in m.iterrows():
        ax.plot(r["resolved"], r["a"], "o", ms=11.0,
                color=BIOME_COLOUR[r["biome"]], mec="white", mew=1.5)
    # least-squares line
    coef = np.polyfit(m["resolved"], m["a"], 1)
    xs = np.linspace(m["resolved"].min() - 0.01,
                     m["resolved"].max() + 0.01, 50)
    ax.plot(xs, np.polyval(coef, xs), "-", color="#5F5E5A", lw=2.0,
            label="least-squares fit")
    r2 = np.corrcoef(m["resolved"], m["a"])[0, 1] ** 2
    ax.text(0.97, 0.94, f"$R^{2}={r2:.2f}$", transform=ax.transAxes,
            fontsize=FS_ANNOT, ha="right", va="top")
    ax.set_xlabel("resolved fraction of the sill,\n"
                  r"$\gamma(h_{\mathrm{trunc}})\,/\,\sigma^{2}$")
    ax.set_ylabel("mean |bias| of the truncated fit (%)")
    ax.set_ylim(0, max(m["a"]) * 1.18)
    ax.set_xlim(m["resolved"].min() - 0.03,
                m["resolved"].max() + 0.02)
    ax.set_title("(a) A truncated fit is biased by exactly\n"
                 "the sill it cannot see", loc="left")
    ax.legend(loc="lower left", fontsize=FS_LEGEND)
    _style(ax)
    # label the two extremes only, to keep the panel readable
    lo_i = m["resolved"].idxmin()
    hi_i = m["resolved"].idxmax()
    for i in (lo_i, hi_i):
        ax.annotate(m.loc[i, "window"].split(" ")[0],
                    xy=(m.loc[i, "resolved"], m.loc[i, "a"]),
                    xytext=(6, 7), textcoords="offset points",
                    fontsize=FS_ANNOT - 2.0, color="#333333")

    # ---- (b) resolved fraction vs fitted range -------------------------
    ax = axs[1]
    for _, r in m.iterrows():
        ax.plot(r["a1_m"], r["resolved"], "o", ms=11.0,
                color=BIOME_COLOUR[r["biome"]], mec="white", mew=1.5)
    # the fit is in 1/a
    inv = 1.0 / m["a1_m"].to_numpy()
    coef2 = np.polyfit(inv, m["resolved"].to_numpy(), 1)
    xs = np.linspace(m["a1_m"].min() * 0.9, m["a1_m"].max() * 1.02, 80)
    ax.plot(xs, np.polyval(coef2, 1.0 / xs), "-", color="#5F5E5A", lw=2.0,
            label=r"fit in $1/a$")
    r2b = np.corrcoef(inv, m["resolved"])[0, 1] ** 2
    ax.text(0.97, 0.94, f"$R^{2}={r2b:.2f}$", transform=ax.transAxes,
            fontsize=FS_ANNOT, ha="right", va="top")
    ax.set_xlabel("fitted range $a$ (m)")
    ax.set_ylabel("resolved fraction of the sill")
    ax.set_xlim(m["a1_m"].min() * 0.80, m["a1_m"].max() * 1.06)
    ax.set_ylim(m["resolved"].min() - 0.035,
                m["resolved"].max() + 0.035)
    ax.set_title("(b) The further the fitted range reaches,\n"
                 "the less the truncation loses", loc="left")
    ax.legend(loc="lower right", fontsize=FS_LEGEND,
              bbox_to_anchor=(0.98, 0.10))
    _style(ax)
    for i in (lo_i, hi_i):
        ax.annotate(m.loc[i, "window"].split(" ")[0],
                    xy=(m.loc[i, "a1_m"], m.loc[i, "resolved"]),
                    xytext=(10, -14) if i == hi_i else (8, 10),
                    textcoords="offset points", ha="right" if i == hi_i else "left",
                    fontsize=FS_ANNOT - 2.0, color="#333333")

    fig.subplots_adjust(left=0.105, right=0.985, top=0.875, bottom=0.155)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=ROOT / "figures_v10" / "fig16.panel-alignment.json",
        overlay_svg=ROOT / "figures_v10" / "fig16.panel-alignment.svg",
    )
    for d in OUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / "fig16_truncation_grading.png")
        fig.savefig(d / "fig16_truncation_grading.pdf")
        fig.savefig(d / "fig16_truncation_grading.svg")
    plt.close(fig)
    print("wrote fig16_truncation_grading.png/.pdf/.svg")
    print(f"  R2(resolved, |bias|) = {r2:.3f};  "
          f"R2(1/a, resolved) = {r2b:.3f}")


def main() -> None:
    apply_style()
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": FS_TICK,
        "axes.titlesize": FS_TITLE,
        "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "savefig.dpi": 400,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    figure15()
    figure16()


if __name__ == "__main__":
    main()
