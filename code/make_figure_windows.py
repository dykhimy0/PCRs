"""make_figure_windows.py -- Figure 4: the eleven reference windows.

The geostatistical comparison is run on eleven 10 m windows of the ESA
WorldCover product, chosen to be different landscapes rather than different
addresses on a map, and spanning seven biomes.  A reader should be able to
judge the claim "seven biomes" without taking it on trust, and to see the
texture that sets each window's variogram range.

Every panel is drawn to one recipe, so that the eleven panels can be compared
at a glance:

    coloured frame   the biome
    title            tile code and site
    sub-label        biome and the window's mean tree cover
    shared scale bar one colour ramp for all eleven panels

Panels are ordered by biome, from temperate through to tropical, which is the
order used in Table 5.1 and everywhere else in the paper.  The twelfth cell of
the 4 x 3 grid carries the biome key and the colour bar, so no panel has to
give up its own space for the legend.

Purpose of the set is stated in the caption: three temperate windows hold the
site fixed and vary the class mix; the cool-temperate, alpine, subtropical and
tropical windows extend the range spectrum at both ends; the second boreal
window turns the boreal biome from a single site into a biome with an
estimable within-biome spread.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from figstyle import apply_style, save, FS_LEGEND, FS_ANNOT  # noqa: E402
from windows_v2 import BIOME_ORDER, WINDOWS, centre, label_of  # noqa: E402

apply_style()

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

BIOME_COLOR = {
    "temperate": "#2F5D8A",
    "cool-temperate": "#3E7C8A",
    "mediterranean": "#B4762A",
    "boreal": "#5B4B8A",
    "alpine": "#7A6A4F",
    "subtropical": "#2E7D5B",
    "tropical": "#A8474B",
}

#: short site names, so the title fits the panel width
SHORT_PLACE = {
    "Black Forest, southwest Germany": "Black Forest, DE",
    "Kraichgau, southwest Germany": "Kraichgau, DE",
    "Swabian Jura, southern Germany": "Swabian Jura, DE",
    "Alentejo montado, Portugal": "Alentejo, PT",
    "Andalusian dehesa, Spain": "Andalusia, ES",
    "Finnish Lapland": "Finnish Lapland",
    "Kainuu, Finland": "Kainuu, FI",
    "Xishuangbanna, Yunnan, China": "Xishuangbanna, CN",
    "Wuyi Shan, Fujian, China": "Wuyi Shan, CN",
    "Bomi, Tibet, China": "Bomi, CN",
    "Xiaoxing'anling, Heilongjiang, China": "Xiaoxing'anling, CN",
}

TREE_CMAP = LinearSegmentedColormap.from_list(
    "treecover", ["#F3EFE6", "#CFE0C4", "#7FAE72", "#2F6B45", "#173F27"])

STEP = 3
DOWNSAMPLE = 4


def load(tile: str, row: int, col: int) -> np.ndarray:
    arr = np.load(DATA / f"worldcover_{tile}_{row}_{col}_3072.npy")
    n = arr.shape[0] // STEP
    a = arr[: n * STEP, : n * STEP]
    tree = (a == 10).astype(np.float32)
    frac = tree.reshape(n, STEP, n, STEP).mean(axis=(1, 3))
    return frac[::DOWNSAMPLE, ::DOWNSAMPLE]


def main() -> None:
    order = sorted(range(len(WINDOWS)),
                   key=lambda i: (BIOME_ORDER.index(WINDOWS[i][3]), i))
    ncol, nrow = 4, 3
    fig = plt.figure(figsize=(7.5, 6.35))

    gs = fig.add_gridspec(nrow, ncol, left=0.015, right=0.985, top=0.965,
                          bottom=0.035, hspace=0.30, wspace=0.075)

    covers = []
    for k, i in enumerate(order):
        tile, row, col, biome, _g, place = WINDOWS[i]
        r, c = divmod(k, ncol)
        ax = fig.add_subplot(gs[r, c])
        frac = load(tile, row, col)
        covers.append((biome, float(frac.mean())))
        ax.imshow(frac, cmap=TREE_CMAP, vmin=0, vmax=1, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_color(BIOME_COLOR[biome])
            s.set_linewidth(1.5)
        ax.set_title(f"{tile}  \u00b7  {SHORT_PLACE[place]}", fontsize=8.4, pad=3.2)
        ax.text(0.5, -0.038, f"{biome}  \u00b7  {frac.mean():.0%} tree",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7.8, color=BIOME_COLOR[biome])

    # ---- twelfth cell: key + scale -------------------------------------
    ax = fig.add_subplot(gs[2, 3])
    ax.axis("off")
    present = [b for b in BIOME_ORDER if any(b == bb for bb, _ in covers)]
    ax.text(0.0, 1.02, "biome (windows)", fontsize=FS_LEGEND, fontweight="bold",
            va="top", ha="left", transform=ax.transAxes)
    for j, b in enumerate(present):
        y = 0.925 - j * 0.088
        n = sum(1 for bb, _ in covers if bb == b)
        ax.add_patch(plt.Rectangle((0.0, y - 0.022), 0.11, 0.052,
                                   facecolor=BIOME_COLOR[b], edgecolor="none",
                                   transform=ax.transAxes, clip_on=False))
        ax.text(0.17, y, f"{b} ({n})", fontsize=7.6, va="center", ha="left",
                color="#1A1A1A", transform=ax.transAxes)

    # the colour bar lives inside that same cell, under the key
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, pos.y0 + 0.028, pos.width * 0.94, 0.0125])
    sm = plt.cm.ScalarMappable(cmap=TREE_CMAP, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_ticks([0, 0.5, 1]); cb.set_ticklabels(["0", "0.5", "1"])
    cb.ax.tick_params(labelsize=7.4, length=2.0, pad=1.2)
    cb.outline.set_linewidth(0.6)
    cax.set_title("tree-cover fraction (30 m)", fontsize=7.4, pad=3, loc="left")

    save(fig, "fig13_reference_windows.png")
    for b, cov in covers:
        print(f"  {b:15s} {cov:.3f}")


if __name__ == "__main__":
    main()
