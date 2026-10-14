# -*- coding: utf-8 -*-
"""
make_figure_metric.py -- the two figures the metric rebuild makes possible.

Neither figure could be drawn on the geographic grid, because both are about
the difference between a distance and a count of cells.

figM1_range_ambiguity
    A fitted range of `a` pixels is a ground distance only once a direction is
    chosen: a 3 x 3 block of 1/12000-degree cells is 27.64 m tall at every
    latitude but 27.64 cos(lat) m wide.  Each window therefore gets a bar
    spanning the two readings, with the metric range marked inside it.  The bar
    is the ambiguity; its length is 1/cos(lat) and reaches 2.5x at 66.7 N.

figM2_anisotropy
    East-west against north-south empirical variograms, in metres.  On the
    metric grid a cell is 10 m in both directions, so any difference between
    the two curves is a property of the landscape.  Panel (b) collects the
    fitted range ratio for all eleven windows.

House style (user-specified, applies to every figure in the paper):
    Times New Roman throughout; large panels; a panel aspect close to square
    rather than a letterbox; generous type so the figure is legible after it is
    scaled into the journal column.

Both panels are read against the same question: is the omnidirectional
variogram the manuscript fits an adequate summary of these fields?
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

import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402

from figstyle import apply_style                     # noqa: E402

# mandatory render-time panel-alignment gate for the multi-panel layout
sys.path.insert(0, str(Path("C:/Users/ASUS/.workbuddy/skills/nature-figure/scripts")))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

RES_M = ROOT / "results_metric"
IMPACT = ROOT / "_metric_vs_geographic_impact.csv"
OUT = ROOT / "figures_v10"

# --- house type sizes, matched to Figure 13 -------------------------------
# Same scale as make_figure13.py / make_figure15_16.py: one type size for the
# whole figure set, sized so it survives being scaled into the journal column.
FS_BASE = 20.0
FS_TICK = 20.0
FS_LABEL = 23.0
FS_TITLE = 22.0
FS_LEGEND = 19.0
FS_ANNOT = 19.0

PANEL_W = 11.6         # inches; wider than the print width so nothing is lost
                       # when the panel is scaled down into the column

BIOME_COLOUR = {
    "boreal": "#3d6b8c", "cool-temperate": "#5b8c6b",
    "temperate": "#7ba05b", "mediterranean": "#c9a227",
    "alpine": "#8c7b9c", "subtropical": "#c2703d", "tropical": "#b5453a",
}


def _style_axes(ax) -> None:
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(FS_TICK)
    ax.xaxis.label.set_fontsize(FS_LABEL)
    ax.yaxis.label.set_fontsize(FS_LABEL)


def _order(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values("lat", ascending=False).reset_index(drop=True)


# ==========================================================================
# figM1 -- the range is not a distance
# ==========================================================================
def fig_range_ambiguity(df: pd.DataFrame) -> None:
    d = _order(df)
    fig, ax = plt.subplots(figsize=(PANEL_W, 6.4))

    xs = np.arange(len(d))
    for i, r in d.iterrows():
        lo, hi = r["range_geo_if_ew_m"], r["range_geo_if_ns_m"]
        ax.plot([i, i], [lo, hi], color=BIOME_COLOUR[r["biome"]],
                lw=11.0, solid_capstyle="butt", alpha=0.45, zorder=2)
        ax.plot(i, r["range_met_m"], marker="o", ms=11.0, color="black",
                mfc="white", mew=2.2, zorder=3)

    ax.set_yscale("log")
    ax.set_ylim(90, 7000)
    ax.set_yticks([100, 300, 1000, 3000])
    ax.get_yaxis().set_major_formatter(
        plt.FuncFormatter(lambda v, _p: f"{v:.0f}"))
    ax.set_xticks(xs)
    ax.set_xticklabels(
        [f"{r['window'].split(' ')[0]}\n{r['biome']}"
         for _, r in d.iterrows()], fontsize=FS_TICK - 4.0, rotation=38,
        ha="right", rotation_mode="anchor")
    ax.set_ylabel("fitted range (m)", fontsize=FS_LABEL)
    ax.grid(axis="y", lw=0.6, alpha=0.30)
    ax.set_axisbelow(True)

    handles = [
        Line2D([], [], color="#555555", lw=11.0, alpha=0.45,
               solid_capstyle="butt",
               label="geographic grid: one range in pixels,\nread as a ground distance"),
        Line2D([], [], color="black", marker="o", ms=11.0, mfc="white",
               mew=2.2, lw=0, label="metric grid: the range, in metres"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=FS_LEGEND)

    ax.set_title("A range in pixels is a distance only after a direction "
                 "is chosen", fontsize=FS_TITLE, loc="left", pad=10)

    worst = d.loc[d["range_geo_ambiguity"].idxmax()]
    iw = int(d["range_geo_ambiguity"].idxmax())
    ax.annotate(
        f"{worst['range_geo_px']:.0f} px = "
        f"{worst['range_geo_if_ew_m']:.0f} m east\u2013west or "
        f"{worst['range_geo_if_ns_m']:.0f} m north\u2013south\n"
        f"({worst['range_geo_ambiguity']:.2f}$\\times$); the metric fit "
        f"gives {worst['range_met_m']:.0f} m",
        xy=(iw, worst["range_met_m"]),
        xytext=(0.30, 0.20), textcoords="axes fraction",
        fontsize=FS_ANNOT - 3.0, ha="left", va="top",
        arrowprops=dict(arrowstyle="->", lw=1.2, color="#333333",
                        connectionstyle="arc3,rad=-0.2"))

    fig.tight_layout()
    for d_ in (OUT, ROOT / "figures"):
        d_.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(d_ / f"figM1_range_ambiguity.{ext}",
                        facecolor="white", edgecolor="none",
                        transparent=False, dpi=400)
    plt.close(fig)
    print("wrote figM1_range_ambiguity.png/.pdf/.svg")


# ==========================================================================
# figM2 -- anisotropy
# ==========================================================================
def fig_anisotropy(df: pd.DataFrame) -> None:
    d = _order(df)
    cur = pd.read_csv(RES_M / "metric_anisotropy_curves.csv")

    # the three windows furthest from isotropic, covering both directions
    dev = (d["aniso_ratio_ew_over_ns"] - 1.0).abs()
    high = d["aniso_ratio_ew_over_ns"] > 1.0
    picks = [d.loc[dev[high].idxmax(), "window"]] if high.any() else []
    seen = set(picks)
    for i in dev.sort_values(ascending=False).index:
        if len(picks) == 3:
            break
        if i not in seen:
            picks.append(d.loc[i, "window"])
            seen.add(i)

    # two stacking panels of near-square aspect, one above the other
    fig, axes = plt.subplots(2, 1, figsize=(PANEL_W, 11.0),
                             gridspec_kw={"height_ratios": [1.35, 1.0],
                                          "hspace": 0.34})

    # ---- (a) directional variograms ---------------------------------
    ax = axes[0]
    for k, lab in enumerate(picks):
        s = cur[cur["window"] == lab]
        sig2 = float(s["sigma2"].iloc[0])
        off = k * 1.02
        ax.plot(s["lag_m"], s["gamma_ew"] / sig2 + off, color="#b5453a",
                lw=2.4, label="east–west" if k == 0 else None)
        ax.plot(s["lag_m"], s["gamma_ns"] / sig2 + off, color="#3d6b8c",
                lw=2.4, ls="--", label="north–south" if k == 0 else None)
        rt = float(d.loc[d["window"] == lab, "aniso_ratio_ew_over_ns"].iloc[0])
        ax.annotate(f"{lab.split(' ')[0]}    ratio {rt:.2f}",
                    xy=(2980, off + 0.06),
                    fontsize=FS_ANNOT, ha="right", va="bottom")
        if k < len(picks) - 1:
            ax.axhline(off + 0.93, lw=0.8, color="#bbbbbb")

    ax.set_xlabel("separation (m)", fontsize=FS_LABEL)
    ax.set_ylabel(r"$\gamma(h)\,/\,\sigma^{2}$  (offset by window)",
                  fontsize=FS_LABEL)
    ax.set_ylim(-0.08, 3.10)
    ax.set_yticks([0, 1.0])
    ax.set_xlim(0, 3000)
    ax.grid(lw=0.6, alpha=0.30)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=FS_LEGEND,
              handlelength=2.4)
    ax.set_title("(a) The three windows furthest from isotropic",
                 fontsize=FS_TITLE, loc="left", pad=10)
    _style_axes(ax)

    # ---- (b) the range ratio for every window -----------------------
    ax = axes[1]
    rt = d["aniso_ratio_ew_over_ns"].values
    cols = [BIOME_COLOUR[b] for b in d["biome"]]
    ax.barh(np.arange(len(d)), rt, color=cols, height=0.72, alpha=0.88)
    ax.axvline(1.0, color="black", lw=1.4)
    ax.axvspan(0.8, 1.25, color="#999999", alpha=0.14, zorder=0)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels([f"{r['window'].split(' ')[0]}  {r['biome']}  "
                        f"{r['lat']:.0f}$^\\circ$"
                        for _, r in d.iterrows()], fontsize=FS_TICK - 1.0)
    ax.set_xlim(0.35, 1.45)
    ax.set_xlabel("range east–west $/$ range north–south", fontsize=FS_LABEL)
    ax.grid(axis="x", lw=0.6, alpha=0.30)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    ax.set_title("(b) All eleven windows", fontsize=FS_TITLE, loc="left", pad=10)
    _style_axes(ax)

    n_out = int(((rt < 0.8) | (rt > 1.25)).sum())
    ax.annotate(f"{n_out} of {len(d)} outside $\\pm$25 %",
                xy=(0.98, 0.05), xycoords="axes fraction",
                fontsize=FS_ANNOT, ha="right")

    fig.tight_layout()
    # measure the rendered panel rectangles before export
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=OUT / "figM2.panel-alignment.json",
        overlay_svg=OUT / "figM2.panel-alignment.svg",
    )
    for d_ in (OUT, ROOT / "figures"):
        d_.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(d_ / f"figM2_anisotropy.{ext}",
                        facecolor="white", edgecolor="none",
                        transparent=False, dpi=400)
    plt.close(fig)
    print("wrote figM2_anisotropy.png/.pdf/.svg")


def main() -> None:
    apply_style()
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": FS_BASE,
        "axes.titlesize": FS_TITLE,
        "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "savefig.dpi": 400,
        # keep text editable in the vector exports
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    df = pd.read_csv(IMPACT)
    fig_range_ambiguity(df)
    fig_anisotropy(df)
    print(f"\n{len(df)} windows; range ambiguity "
          f"{df['range_geo_ambiguity'].min():.2f}x to "
          f"{df['range_geo_ambiguity'].max():.2f}x")


if __name__ == "__main__":
    main()
