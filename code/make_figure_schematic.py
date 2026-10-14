"""
make_figure_schematic.py -- Figure 2: the representational gap and the three
information levels.

Panel (a) draws the aggregation step: fine units carrying marginals are pooled
into a block whose law needs the (unstored) joint.  Panel (b) is the metadata
ladder, showing what each stored level lets one answer for mean, variance and
exceedance (identified / bounded / not identified).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

from figstyle import apply_style

apply_style()

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

C_MARG = "#2F5D8A"
C_JOINT = "#A8474B"
GREEN, AMBER, RED = "#3E7C59", "#D1823C", "#A32D2D"


def small_marginal(ax, x, y, w, p, h=0.5):
    """A mini two-part bar hinting at a marginal law."""
    ax.add_patch(Rectangle((x, y), w * p, h, facecolor=C_MARG, edgecolor="none"))
    ax.add_patch(Rectangle((x + w * p, y), w * (1 - p), h,
                           facecolor="#B5D4F4", edgecolor="none"))


def main():
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.4, 4.1),
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    # ---- (a) aggregation ----
    ax = axa
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title("(a) what is stored vs what is needed", fontsize=12.5)

    for r in range(2):
        for c in range(2):
            x, y = 0.7 + c * 1.9, 5.6 + r * 1.9
            ax.add_patch(FancyBboxPatch((x, y), 1.5, 1.5,
                                        boxstyle="round,pad=0.04,rounding_size=0.12",
                                        linewidth=1.0, edgecolor=C_MARG, facecolor="#EDF3FA"))
            small_marginal(ax, x + 0.32, y + 0.35, 0.85, 0.35 + 0.05 * ((r + c) % 3))
            ax.text(x + 0.75, y + 1.16, r"$\Lambda_j$", ha="center", va="bottom",
                    fontsize=11, color=C_MARG)
    ax.text(1.7, 4.7, "fine units carry marginals", ha="center", fontsize=11,
            color=C_MARG)

    ax.add_patch(FancyArrowPatch((4.1, 6.35), (6.1, 6.35), arrowstyle="-|>",
                                 mutation_scale=22, linewidth=1.6, color="#1A1A1A"))
    ax.text(5.1, 6.75, r"$g_w$", ha="center", fontsize=12)

    ax.add_patch(FancyBboxPatch((6.3, 5.15), 3.1, 2.4,
                                boxstyle="round,pad=0.06,rounding_size=0.14",
                                linewidth=1.2, edgecolor=C_JOINT, facecolor="#FBEFF0"))
    ax.text(7.85, 6.95, r"$\Lambda_C = (g_w)_\#\,\Gamma$", ha="center",
            va="top", fontsize=12.5, color=C_JOINT)
    ax.text(7.85, 5.95, "needs the joint law", ha="center", fontsize=11, color=C_JOINT)
    ax.text(7.85, 5.45, "(not stored)", ha="center", fontsize=11,
            style="italic", color=C_JOINT)

    # ---- (b) ladder ----
    ax = axb
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title("(b) the metadata ladder", fontsize=12.5)

    cols = ["stored", "mean", "variance", "exceedance"]
    xs = [0.6, 3.4, 5.9, 8.4]
    ax.text(xs[0], 9.2, "level", fontsize=11, fontweight="bold", color="#1A1A1A")
    for c, name in zip(cols[1:], ["mean", "variance", "exceedance"]):
        ax.text(xs[cols.index(name)] + 0.6, 9.2, name, ha="center",
                fontsize=11, fontweight="bold", color="#1A1A1A")

    rows = [
        ("A", "marginals", GREEN, AMBER, RED),
        ("B", "+ pair covariance", GREEN, GREEN, AMBER),
        ("C", "full joint", GREEN, GREEN, GREEN),
    ]
    MARK = {GREEN: "\u25cf", AMBER: "\u25cb", RED: "\u00d7"}
    for i, (lvl, store, mstat, vstat, estat) in enumerate(rows):
        y = 8.0 - i * 2.0
        ax.add_patch(FancyBboxPatch((0.5, y - 0.55), 2.4, 1.5,
                                    boxstyle="round,pad=0.04,rounding_size=0.1",
                                    linewidth=0.8, edgecolor="#5F5E5A", facecolor="#F4F3EF"))
        ax.text(1.7, y + 0.2, lvl, ha="center", fontsize=13,
                fontweight="bold", color="#1A1A1A")
        ax.text(1.7, y - 0.35, store, ha="center", fontsize=10, color="#444441")
        for stat, cname in zip((mstat, vstat, estat), ("mean", "variance", "exceedance")):
            cx = xs[cols.index(cname)] + 0.6
            ax.text(cx, y + 0.2, MARK[stat], ha="center", va="center",
                    fontsize=16, color=stat)

    # legend for the dots
    for yy, (ch, col, txt) in zip((1.7, 1.1, 0.5), [
            ("\u25cf", GREEN, "identified"),
            ("\u25cb", AMBER, "bounded / partial"),
            ("\u00d7", RED, "not identified")]):
        ax.text(2.0, yy, ch, fontsize=13, color=col, va="center")
        ax.text(2.6, yy, txt, fontsize=10.5, color="#1A1A1A", va="center")

    fig.tight_layout()
    fig.savefig(FIG / "fig0_schematic.png")
    plt.close(fig)
    print("wrote fig0_schematic.png")


if __name__ == "__main__":
    main()
