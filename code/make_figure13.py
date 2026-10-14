# -*- coding: utf-8 -*-
"""
make_figure13.py -- asymptotic sharpness of the PSD relaxation (Figure 13).

Redrawn to the paper's house style and to the corrected terminology.

Terminology
-----------
Panel (b) plots gap / sigma^2, an **absolute gap normalised by the single-unit
marginal variance**.  It is NOT a relative error: when the relaxation returns a
zero gap while the exact lower endpoint is positive, a relative error against
the exact endpoint is undefined.  The axis label, the panel title and the legend
therefore say "gap / sigma^2" and "normalised gap", never "relative error".

Style
-----
Times New Roman throughout (mathtext in STIX so the symbols match the serif
face); panels laid out large and nearly square; type sized so that it stays
legible after scaling into the journal column.

Data: results/psd_sweep_m.csv (panel a, b), closed forms from exp_psd.py (panel c).
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

from exp_psd import gap_closed, gap_envelope, rel_gap_envelope  # noqa: E402
from figstyle import apply_style                      # noqa: E402

# mandatory render-time panel-alignment gate for multi-panel layouts
sys.path.insert(0, str(Path("C:/Users/ASUS/.workbuddy/skills/nature-figure/scripts")))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402

RES = ROOT / "results"
OUT_DIRS = (ROOT / "figures_v10", ROOT / "figures", ROOT / "figures_v2")

# --- house type sizes, sized for a large canvas --------------------------
# Raised deliberately: the panels are wide enough that the type can be large
# without any label collision, and the figure must stay readable when it is
# scaled down into the journal column.
FS_TICK = 20.0
FS_LABEL = 23.0
FS_TITLE = 22.0
FS_LEGEND = 19.0
FS_ANNOT = 19.0

C_ENV = "#5F5E5A"
C_05 = "#2F5D8A"
C_03 = "#B03A5B"
C_IRR = "#3E7C59"
C_REF = "#B4B2A9"


def _style(ax) -> None:
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(FS_TICK)
    ax.xaxis.label.set_fontsize(FS_LABEL)
    ax.yaxis.label.set_fontsize(FS_LABEL)
    ax.title.set_fontsize(FS_TITLE)


def figure13() -> None:
    df = pd.read_csv(RES / "psd_sweep_m.csv")

    # Two rows: (a) and (b) share the m axis on top, (c) spans the bottom.
    # A single 3-wide row made every panel narrow and the type cramped; the
    # 2 + 1 layout gives each panel roughly twice the plot area at the same
    # type size.  White background throughout.
    fig = plt.figure(figsize=(14.0, 13.0), facecolor="white")
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[1.0, 0.92],
        left=0.075, right=0.985, top=0.945, bottom=0.062,
        wspace=0.24, hspace=0.30,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    axs = [ax_a, ax_b, ax_c]

    # ---- (a) gap vs m, log-log, with the 1/(4 m^2) envelope ------------
    ax = axs[0]
    ms = np.arange(2, 65)
    # the envelope first, so the data sit on top of it
    ax.plot(ms, gap_envelope(ms), "-", color=C_ENV, lw=2.8,
            label="envelope $1/(4m^{2})$", zorder=1)
    # the gap oscillates between 0 and the envelope as frac(m p) cycles, so the
    # closed form is drawn as markers only; joining them draws a spurious sawtooth
    for p, col, mk, lab in ((0.5, C_05, "o", "$p=0.5$"),
                            (0.3, C_03, "s", "$p=0.3$"),
                            (1 / np.pi, C_IRR, "^", r"$p=1/\pi$")):
        g = df[df["p"].round(12) == round(p, 12)].sort_values("m")
        gg = g[g["gap"] > 0]
        ax.plot(gg["m"], gg["gap"], linestyle="none", marker=mk, ms=6.0,
                mfc="white", mew=1.4, color=col, label=lab, zorder=3)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.set_ylim(1e-9, 3.0)
    ax.set_xlabel("pixels aggregated, $m$")
    ax.set_ylabel("joint-mixability gap")
    ax.set_title("(a) The gap is bounded by $1/(4m^{2})$", loc="left", pad=12)
    ax.legend(loc="lower left", fontsize=FS_LEGEND, handlelength=2.0,
              labelspacing=0.55)
    ax.grid(lw=0.5, alpha=0.25)
    ax.set_axisbelow(True)
    _style(ax)

    # ---- (b) normalised gap vs m ---------------------------------------
    # NOTE: this is gap / sigma^2, an absolute gap normalised by the
    # single-unit marginal variance -- explicitly not a relative error.
    ax = axs[1]
    ax.plot(ms, rel_gap_envelope(0.5, ms), "-", color=C_05, lw=2.4,
            label="$1/m^{2}$   ($p=0.5$)", zorder=2)
    ax.plot(ms, rel_gap_envelope(0.3, ms), "-", color=C_03, lw=2.1,
            label="$1/(4m^{2}p(1-p))$   ($p=0.3$)", zorder=2)
    for p, col in ((0.5, C_05), (0.3, C_03)):
        g = df[(df["p"].round(12) == round(p, 12))].sort_values("m")
        ax.plot(g["m"], g["rel_gap"].clip(lower=1e-9), "o", ms=6.0,
                color=col, alpha=0.90, zorder=3)
    ax.axhline(0.01, color=C_REF, lw=1.3, ls=":")
    ax.text(2.18, 0.0128, "1 %", fontsize=FS_ANNOT - 0.5, color=C_ENV)
    ax.axvline(10, color=C_REF, lw=1.1, ls=":")
    # place the m = 10 marker inside the axes: below the title row and clear
    # of the legend in the lower left
    ax.text(11.0, 0.012, "$m=10$", fontsize=FS_ANNOT - 2.0, color=C_ENV,
            rotation=90, va="top", ha="left", rotation_mode="anchor")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.set_ylim(1e-5, 0.5)
    ax.set_xlabel("pixels aggregated, $m$")
    # corrected label: the axis is a gap normalised by the marginal variance
    ax.set_ylabel("normalised gap,  gap $/\\,\\sigma^{2}$")
    ax.set_title("(b) The normalised gap falls as $m^{-2}$", loc="left", pad=12)
    ax.legend(loc="lower left", fontsize=FS_LEGEND, handlelength=2.0,
              labelspacing=0.55)
    ax.grid(lw=0.5, alpha=0.25)
    ax.set_axisbelow(True)
    _style(ax)

    # ---- (c) gap vs p at fixed m -- the two regimes --------------------
    ax = axs[2]
    ps = np.linspace(0.005, 0.995, 1600)
    cols = {4: "#1F4E79", 8: "#2E8B7A", 16: "#B03A5B"}
    # shade the region where m min(p, 1-p) > 1, i.e. where the PSD corner binds
    lim16 = 16 * np.minimum(ps, 1 - ps)
    span = ps[lim16 > 1]
    if span.size:
        ax.axvspan(span.min(), span.max(), color="#EFEBE0", zorder=0)
    for m in (4, 8, 16):
        g = np.array([gap_closed(float(p), m) for p in ps])
        ax.plot(ps, g, "-", lw=2.4, color=cols[m], label=f"$m={m}$", zorder=3)
        ax.axhline(gap_envelope(m), color=cols[m], lw=1.2, ls=":", zorder=2)
    # the two regime boundaries, m min(p, 1-p) = 1
    for m in (4, 8, 16):
        pk = 1.0 / m
        for xb in (pk, 1 - pk):
            if 0 < xb < 1:
                ax.plot([xb], [0], marker="v", ms=7.0, color=cols[m],
                        clip_on=False, zorder=4)
    ax.set_xlabel("class proportion $p$")
    ax.set_ylabel("joint-mixability gap")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.0008, 0.0245)
    ax.set_title("(c) Zero until $m\\min(p,1-p)>1$ (marked $\\nabla$), then bounded by $1/(4m^{2})$",
                 loc="left", pad=12)
    # sit below the legend row, which occupies the top of the axes
    ax.annotate("shaded: the PSD corner binds",
                xy=(0.5, 0.0192), xycoords=("data", "data"),
                fontsize=FS_ANNOT - 2.0, color="#6B6B6B",
                ha="center", va="center")
    ax.legend(loc="upper left", ncol=3, columnspacing=1.1,
              handlelength=1.6, fontsize=FS_LEGEND)
    ax.grid(lw=0.5, alpha=0.25)
    ax.set_axisbelow(True)
    _style(ax)

    # mandatory: measure the final rendered panel rectangles before export
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=ROOT / "figures_v10" / "fig13.panel-alignment.json",
        overlay_svg=ROOT / "figures_v10" / "fig13.panel-alignment.svg",
    )

    for d in OUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf", "svg"):
            # white background, opaque: the panels are read on the printed page
            fig.savefig(d / f"fig13_psd_asymptotic.{ext}",
                        facecolor="white", edgecolor="none",
                        transparent=False, dpi=400)
    plt.close(fig)
    print("wrote fig13_psd_asymptotic.png/.pdf/.svg ->",
          ", ".join(d.name for d in OUT_DIRS))


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
        # keep text editable in both vector formats
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    figure13()


if __name__ == "__main__":
    main()
