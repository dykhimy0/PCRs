"""make_figure_sets.py -- Figure 3: the geometry of the identified set.

Section 4 proves three things that are easy to state and hard to picture:
the pairwise (Frechet-Hoeffding) relaxation can leave the feasible set, the
positive-semidefinite relaxation never does, and which of the two the exact
program agrees with is decided by a single inequality in the block size and
the class proportion.

Panel (a) draws the three intervals for one concrete case, so that the nesting
exact  subset  PSD  subset  FH is visible and the negative FH endpoint is seen
rather than asserted.  Panel (b) draws the regime boundary itself for every
block size and proportion a practitioner might meet.

Revision v10 makes three changes to the v9 figure:

* the two floating text blocks in panel (a) ("exact = PSD ...", "negative: not
  a variance") are gone; they covered the bars, so their content now lives in
  a four-entry legend above the bars, and m and p moved into the title;
* the same label was also numerically wrong.  At m = 4 and p = 0.30 the regime
  predicate gives m*min(p, 1-p) = 1.2 > 1, i.e. a strictly POSITIVE gap, and
  indeed exact = 0.0100 while PSD = 0.0000.  The v9 text asserted equality;
  the v10 legend states the three endpoints as computed;
* in panel (b) the "relaxation exact" label sat at m = 30, p = 0.055, which
  the panel's own boundary puts in the NON-exact region (1/p = 18.2 < 30).
  It is moved to m = 6.5, p = 0.075, where m*p = 0.49 <= 1.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from figstyle import apply_style, save  # noqa: E402

apply_style()

ROOT = Path(__file__).resolve().parent.parent

C_EXACT = "#2F5D8A"
C_PSD = "#7B5EA7"
C_FH = "#B4762A"
C_BAD = "#9A3A3A"
C_OK = "#DCE9F5"
C_GAP = "#F6E3D2"


def sharp_min(m: int, p: float) -> float:
    """Lower endpoint of the sharp interval (Proposition 5)."""
    pi = (m * p) % 1.0
    return pi * (1.0 - pi) / m ** 2


def psd_min(m: int, p: float) -> float:
    """Lower endpoint of the PSD relaxation (Proposition 6)."""
    s2 = p * (1.0 - p)
    rho_fh = -min(p, 1.0 - p) / max(p, 1.0 - p)
    rho = max(rho_fh, -1.0 / (m - 1))
    return s2 * (1.0 / m + (m - 1.0) / m * rho)


def fh_min(m: int, p: float) -> float:
    """Lower endpoint of the pairwise Frechet-Hoeffding relaxation."""
    s2 = p * (1.0 - p)
    rho_fh = -min(p, 1.0 - p) / max(p, 1.0 - p)
    return s2 * (1.0 / m + (m - 1.0) / m * rho_fh)


def main() -> None:
    m0, p0 = 4, 0.30
    up = p0 * (1 - p0)
    lo_exact, lo_psd, lo_fh = sharp_min(m0, p0), psd_min(m0, p0), fh_min(m0, p0)

    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.4),
                            gridspec_kw={"width_ratios": [1.0, 1.12]})

    # ---------------- (a) the three intervals --------------------------
    ax = axs[0]
    rows = [
        ("Fr\u00e9chet\u2013Hoeffding\nrelaxation", lo_fh, C_FH),
        ("PSD relaxation", lo_psd, C_PSD),
        ("exact linear program", lo_exact, C_EXACT),
    ]
    ypos = np.arange(len(rows))[::-1]
    ax.axvline(up, color="#8A8A83", lw=1.0, ls="--", zorder=1)
    ax.text(up, -0.44, r"sill $=\sigma^2$", ha="center", va="bottom",
            fontsize=8.6, color="#6E6E68")

    for y, (lab, lo, col) in zip(ypos, rows):
        ax.plot([max(lo, 0.0), up], [y, y], "-", lw=7.0, color=col,
                solid_capstyle="butt", alpha=0.90, zorder=3)
        ax.plot([lo, lo], [y - 0.20, y + 0.20], "-", lw=1.4, color=col,
                zorder=4)
        if lo < 0:
            ax.plot([lo, 0.0], [y, y], "-", lw=2.0, color=C_BAD, zorder=4)
            ax.plot([lo], [y], "x", ms=7, mew=1.8, color=C_BAD, zorder=5)
        # the lower endpoint is the quantity that carries the message, so it
        # is printed under the left end of its own bar rather than at the sill
        ax.text(lo, y - 0.30, f"{lo:.3f}", ha="center", va="top",
                fontsize=8.8, color=col)

    ax.set_yticks(ypos)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.0)
    ax.set_xlim(-0.055, up + 0.10)
    ax.set_ylim(-0.85, 3.95)
    ax.set_xlabel(r"lower endpoint of the interval for $\mathrm{Var}(\nu_C)$",
                  fontsize=9.4)
    ax.set_title(r"(a) The three intervals do not coincide"
                 "\n" r"($m=4$, $p=0.30$)", fontsize=10.4)

    handles = [
        Line2D([], [], color=C_EXACT, lw=5.0, label="exact linear program"),
        Line2D([], [], color=C_PSD, lw=5.0, label="PSD relaxation"),
        Line2D([], [], color=C_FH, lw=5.0,
               label="Fr\u00e9chet\u2013Hoeffding relaxation"),
        Line2D([], [], color=C_BAD, lw=1.8, ls=(0, (4, 2)), marker="x",
               ms=6, mew=1.6, label="endpoint below zero: not a variance"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.01),
              ncol=1, fontsize=8.4, frameon=False, handlelength=1.6,
              labelspacing=0.30, borderaxespad=0.0)

    # ---------------- (b) the regime boundary --------------------------
    ax = axs[1]
    ms = np.unique(np.round(np.logspace(np.log10(2), np.log10(64), 260)).astype(int))
    ps = np.linspace(0.005, 0.5, 240)
    M, P = np.meshgrid(ms, ps)
    exact = (M * np.minimum(P, 1 - P)) <= 1.0

    # Bands are keyed to the predicate, not to the runner's taste: the cool
    # fill is where the relaxation is exact, the warm fill where the gap is
    # positive.  v9 passed the two colours in the opposite order, so its two
    # in-plot labels named the wrong regions.
    ax.contourf(M, P, exact.astype(float), levels=[-0.5, 0.5, 1.5],
                colors=[C_GAP, C_OK], zorder=0)
    ax.plot(1.0 / ps, ps, "-", lw=1.9, color="#B4762A", zorder=3)

    ax.set_xscale("log")
    ax.set_xlim(2, 64)
    ax.set_ylim(0, 0.5)
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.minorticks_off()
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xlabel("pixels per block, $m$", fontsize=9.4)
    ax.set_ylabel("within-block proportion, $p$", fontsize=9.4)
    ax.set_title(r"(b) The boundary is $m\min(p,1-p)=1$", fontsize=10.4)

    # the region names were in-plot text in v9 and crossed their own boundary;
    # they are a legend now, which is why the two labels sit in one box
    handles_b = [
        Line2D([], [], color="#B4762A", lw=1.9,
               label=r"boundary $m\min(p,1-p)=1$"),
        Patch(facecolor=C_OK, edgecolor="none",
              label="pairwise relaxation = exact"),
        Patch(facecolor=C_GAP, edgecolor="none",
              label="positive joint-mixability gap"),
    ]
    ax.legend(handles=handles_b, loc="upper right", fontsize=8.3,
              frameon=False, labelspacing=0.32, handlelength=1.6,
              borderaxespad=0.15)
    ax.text(-0.015, -0.265,
            r"at $m=2$ the pairwise relaxation is exact for every $p$",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.3,
            color="#444441")

    fig.subplots_adjust(left=0.205, right=0.985, top=0.86, bottom=0.255,
                        wspace=0.50)
    save(fig, "fig12_identified_set_geometry.png")
    print(f"  m={m0} p={p0}: FH={lo_fh:.4f}  PSD={lo_psd:.4f}  "
          f"exact={lo_exact:.4f}  upper={up:.4f}  "
          f"gap={lo_exact:.4f} <= 1/(4m^2)={1/(4*m0**2):.4f}")


if __name__ == "__main__":
    main()
