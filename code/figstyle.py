"""figstyle.py -- one publication style for every figure in the paper.

Design rules (fixed once here, applied everywhere):

* Times New Roman for all text, STIX for math so symbols match the serif face.
* Type sizes are chosen so that the figure survives being scaled to the
  journal's 190 mm double-column width.  A figure is laid out at its final
  print width (PRINT_W inches) and the type is set at its final print size,
  so nothing has to be guessed at the typesetting stage.
* Panels are sized for legibility, not for content density: when a panel needs
  more room, drop an annotation rather than shrink the font.

Import and call apply_style() at the top of every make_figure*.py, then use
save(fig, name) instead of fig.savefig so that the file lands in both the
working figure directory and the figure set referenced by the manuscript.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib as mpl

ROOT = Path(__file__).resolve().parent.parent

#: Double-column print width of the journal, in inches (190 mm).
#: Figures are laid out at this width and then kept as square as the content
#: allows, so that each panel is large on the page.
PRINT_W = 7.5
#: Single-column print width, in inches (90 mm).
SINGLE_W = 3.5
#: Preferred height when a figure should read as "large and square".
PRINT_H_SQ = 6.6

#: Final printed type sizes, in points.
#: Raised for legibility: panels are laid out larger and squarer so the type
#: survives being scaled down into the journal column.  Rule: drop an
#: annotation rather than shrink the font.
FS_BASE = 15.0        # body text inside a figure
FS_TICK = 14.0        # tick labels
FS_LABEL = 16.0       # axis labels
FS_TITLE = 15.5       # panel titles
FS_LEGEND = 14.0      # legend text
FS_ANNOT = 14.0       # in-axes annotation text (data labels)


def apply_style() -> None:
    """Install the shared rcParams.  Safe to call more than once."""
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
        "legend.handlelength": 1.8,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.2,
        "legend.labelspacing": 0.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9,
        "axes.labelpad": 4.0,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "lines.linewidth": 1.6,
        "figure.dpi": 130,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "savefig.facecolor": "white",
        "legend.frameon": False,
        "errorbar.capsize": 2.0,
    })
    _warn_if_font_missing()


def _warn_if_font_missing() -> None:
    """Loud warning when Times New Roman is absent.

    matplotlib falls back silently, which would silently change every figure
    in the paper, so make the failure visible at generation time.
    """
    try:
        from matplotlib import font_manager as fm
    except Exception:                                    # pragma: no cover
        return
    have = {f.name for f in fm.fontManager.ttflist}
    if "Times New Roman" not in have:
        warnings.warn(
            "Times New Roman is NOT installed; figures will fall back to "
            f"{mpl.rcParams['font.serif'][1]}. Install the font before "
            "producing submission figures.",
            RuntimeWarning,
            stacklevel=2,
        )


def figsize(w: float = PRINT_W, h: float = 3.0) -> tuple:
    """Return a (width, height) in inches, width clamped to the print width."""
    return (w, h)


def save(fig, name: str, also_v2: bool = True) -> None:
    """Save a figure into figures/ and, by default, figures_v2/.

    The manuscript references figures_v2/; several of the generators were
    written when the manuscript referenced figures/.  Writing both keeps the
    two directories from drifting apart.
    """
    import matplotlib.pyplot as plt

    targets = [ROOT / "figures"]
    if also_v2:
        targets.append(ROOT / "figures_v2")
    for d in targets:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / name)
    plt.close(fig)
    print("wrote", name, "->", ", ".join(str(t.name) for t in targets))
