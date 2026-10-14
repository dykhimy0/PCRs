"""
make_figure_roadmap.py -- Figure 1: research roadmap.

A single-column vertical flowchart of the study, from the motivating question
to the decision-level consequence.  Hand-drawn with matplotlib patches so the
type is Times New Roman and sized for print.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from figstyle import apply_style

apply_style()

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# each stage: (title, body)
STAGES = [
    ("1. Question",
     "Does the information a probabilistic raster stores determine the law of a "
     "coarser unit derived from it?"),
    ("2. Formalisation",
     "Probabilistic-compositional pixel (PCR); aggregate = area-weighted push-forward; "
     "three stored levels A (marginals), B (+pair covariance), C (full joint)."),
    ("3. Theory",
     "Identifiability of mean, variance, exceedance; sharp bounds by linear "
     "programming; Fréchet\u2013Hoeffding and PSD relaxations (Prop. 1\u20136)."),
    ("4. Experiments",
     "Synthetic landscapes with known truth; exact-LP verification; six ESA "
     "WorldCover windows \u00d7 three biomes; geostatistical comparator; decision cost."),
    ("5. Findings",
     "Mean is safe; variance needs pair information; independence understates block "
     "variance 3.5\u201319.8\u00d7; store dependence or store the sharp bounds."),
]

C_STAGE = "#2F5D8A"
C_EDGE = "#5F5E5A"


def main():
    fig, ax = plt.subplots(figsize=(6.6, 8.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 15)
    ax.axis("off")

    n = len(STAGES)
    x0, x1 = 0.6, 9.4
    y_top, y_bot = 13.6, 1.4
    gap = 0.6
    box_h = (y_top - y_bot - (n - 1) * gap) / n

    for i, (title, body) in enumerate(STAGES):
        y1 = y_top - i * (box_h + gap)
        y0 = y1 - box_h
        box = FancyBboxPatch((x0, y0), x1 - x0, box_h,
                             boxstyle="round,pad=0.08,rounding_size=0.15",
                             linewidth=1.0, edgecolor=C_STAGE, facecolor="#EDF3FA")
        ax.add_patch(box)
        ax.text((x0 + x1) / 2, y0 + box_h - 0.28, title,
                ha="center", va="top", fontsize=13, fontweight="bold", color=C_STAGE)
        wrapped = textwrap.fill(body, width=62)
        ax.text((x0 + x1) / 2, y0 + 0.22, wrapped,
                ha="center", va="bottom", fontsize=11, color="#1A1A1A",
                linespacing=1.35)

        if i < n - 1:
            y_from = y0
            y_to = y_from - gap + 0.04
            ax.add_patch(FancyArrowPatch(((x0 + x1) / 2, y_from),
                                         ((x0 + x1) / 2, y_to),
                                         arrowstyle="-|>", mutation_scale=18,
                                         linewidth=1.4, color=C_EDGE))

    fig.tight_layout()
    fig.savefig(FIG / "fig0_roadmap.png")
    plt.close(fig)
    print("wrote fig0_roadmap.png")


if __name__ == "__main__":
    main()
