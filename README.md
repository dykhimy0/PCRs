# Scale-consistency of probabilistic-compositional rasters (PCRs)

Reproducibility repository for the manuscript

> **Scale-consistency of probabilistic-compositional rasters: what a per-pixel
> distribution does and does not determine**

Everything needed to regenerate every number, table and figure in the paper is
in this repository. The repository is self-contained apart from the ESA
WorldCover 2021 v200 tiles, which are public and are re-downloadable (see
[Data](#3-data)).

---

## 1 What the study does

A probabilistic-compositional raster (PCR) stores a probability distribution
over the class fractions inside each fine cell. The paper asks which properties
of a **coarser** unit are determined by that stored information, and answers it
as a partial-identification problem over three levels of stored information:

| Level | Stored | Aggregate mean | Aggregate variance | Exceedance |
|---|---|---|---|---|
| **A** | marginals `Lambda_j` | identified | **not** identified | not identified |
| **B** | marginals + pair covariance | identified | identified exactly | bounded |
| **C** | full joint `Gamma` | identified | identified | identified |

The pivot is **Corollary 2.1**: the aggregate variance needs pair covariance and
nothing more; a per-pixel variance field stores exactly the part that does not
suffice. The size of the resulting ambiguity is the **aggregation-ambiguity
index** (ACI), `(V+ - V-) / V+`, bounded below by `1 - 1/m`.

---

## 2 Repository layout

```
code/                 analysis and figure scripts (see below)
results/              CSV outputs of the theory + synthetic experiments
results_metric/       CSV outputs of the metric-support real-data analysis
figures/              the manuscript figures (PNG, 400 dpi)
docs/                 window metadata and the figure-order reference
environment.yml       pinned package versions
SHA256SUMS            checksums of every tracked file
```

### `code/` — core library

| File | Role |
|---|---|
| `pcp_core.py` | transportation-polytope construction, exact LP (HiGHS via SciPy), FH and PSD relaxations, closed forms |
| `windows_v2.py` | **single source of truth** for the eleven analysis windows (tile, row, col, biome, place, centre latitude) |
| `window_metadata.py` | per-window latitude/longitude, ground dimensions, biome key |
| `metric_fields.py` | azimuthal-equidistant reprojection to a common metric support |
| `figstyle.py` | shared matplotlib house style (Times New Roman, STIX mathtext, type sizes) |

### `code/` — experiments

| Script | Produces | Read by |
|---|---|---|
| `exp_theory.py` | `identification_ladder.csv`, `exceedance_envelope.csv`, `metadata_ladder.csv`, `copula_family.csv` | Fig 5, 6 |
| `exp_spatial.py` | `landscape_sweep.csv`, `coverage_calibration.csv` | Fig 7, 8 |
| `exp_real.py`, `exp_real_v2.py` | `real_*.csv` | Fig 9 |
| `exp_baseline.py` | `baseline_*.csv` | Fig 15 (geographic-grid version) |
| `exp_compositional.py` | `compositional_*.csv` | Fig 14 |
| `exp_estimator.py` | `estimator_shrinkage.csv` | Fig 12 |
| `exp_psd.py` | `psd_sweep_m.csv`, `psd_fh_table.csv`, `pairwise_relaxation_gap.csv` | Fig 13 |
| `exp_metric.py` | all of `results_metric/` | Fig 15, 16, 17, 18 |

### `code/` — figures

`make_figure_roadmap.py` (Fig 1), `make_figure_schematic.py` (Fig 2),
`make_figure_sets.py` (Fig 3), `make_figure_windows.py` (Fig 4),
`make_figures.py` (Fig 5, 6, 8, 9, 10), `make_figure7.py` (Fig 7),
`make_figure8.py` (Fig 14), `make_figure_rho.py` (Fig 11),
`make_figure_estimator_v12.py` (Fig 12), `make_figure13.py` (Fig 13),
`make_figure15_16.py` (Fig 15, 16), `make_figure_metric.py` (Fig 17, 18).

---

## 3 Data

### Synthetic

Fully reproducible from `exp_*.py`. Every generator uses an explicit NumPy
`default_rng` seed; no external input is required.

### ESA WorldCover 2021 v200

The eleven windows are sub-windows of public ESA WorldCover 2021 v200 tiles.
Each window is `3072 x 3072` cells of the native `1/12000` degree grid.

- **Product**: ESA WorldCover 2021 v200, 10 m, 3 x 3 degree tiles, 36000 x 36000 cells.
- **Access**: cloud-optimised GeoTIFF,
  `s3://esa-worldcover/v200/2021/map/ESA_WorldCover_10m_2021_v200_<tile>_Map.tif`
  (anonymous, no credentials). Public mirrors and the AWS open-data registry
  both serve the same objects.
- **Tile naming**: a tile is named after its **south-west** corner and runs
  3 degrees north and east. Row 0 is at the **north** edge, so a window whose
  top-left cell is `(row, col)` has its centre at

  ```
  lat = tile_N + 3 - (row + 1536) / 12000
  lon = tile_E     + (col + 1536) / 12000
  ```

  This convention is implemented in `code/windows_v2.py`.

- **Class of interest**: tree cover (class 10) is mapped to a binary indicator;
  all other classes are pooled. For the multi-class experiment, tree (10),
  crop (40) and pooled-other are projected onto the compositional lattice.

The eleven windows, and the underlying tiles to download:

| Window | Tile | Biome | Site |
|---|---|---|---|
| N48E006 (26k, 26k) | `N48E006` | temperate | Black Forest, Germany |
| N48E006 (28k, 20k) | `N48E006` | temperate | Kraichgau, Germany |
| N48E009 (26k, 26k) | `N48E009` | temperate | Swabian Jura, Germany |
| N39W009 (18k, 18k) | `N39W009` | mediterranean | Alentejo, Portugal |
| N36W006 (15k, 15k) | `N36W006` | mediterranean | Andalusian dehesa, Spain |
| N66E024 (26k, 20k) | `N66E024` | boreal | Finnish Lapland |
| N63E027 (14k, 6k) | `N63E027` | boreal | Kainuu, Finland |
| N21E099 (24k, 12k) | `N21E099` | tropical | Xishuangbanna, China |
| N24E117 (14k, 6k) | `N24E117` | subtropical | Wuyi Shan, China |
| N27E093 (10k, 16k) | `N27E093` | alpine | Bomi, Tibet, China |
| N45E126 (2k, 25k) | `N45E126` | cool-temperate | Xiaoxing'anling, China |

Eight distinct tiles are required.

### Metric-support reprojection

All variogram and block-variance statistics in Sections 5.3–5.6 are computed on
a **common metric support**, because a cell of a geographic grid is not a fixed
ground distance. Each window is reprojected with a per-window
azimuthal-equidistant projection:

```
+proj=aeqd +lat_0=<window centre> +lon_0=<window centre> +datum=WGS84 +units=m
```

with `calculate_default_transform(..., resolution=10)` and
`Resampling.nearest`, recoding to the binary indicator **before** resampling so
that no interpolation mixes class labels. On this grid one cell is 10 m in both
directions, so any east–west / north–south difference is a property of the
landscape rather than of the grid.

---

## 4 How to run

```bash
conda env create -f environment.yml
conda activate pcr

# 1. theory and synthetic experiments
python code/exp_theory.py
python code/exp_spatial.py
python code/exp_psd.py
python code/exp_compositional.py
python code/exp_estimator.py

# 2. real-data experiments (needs the WorldCover tiles in data/)
python code/exp_real_v2.py
python code/exp_metric.py          # metric support: results_metric/

# 3. figures
python code/make_figure13.py
python code/make_figure15_16.py
python code/make_figure_metric.py
```

Every script prints the numbers it produces and writes CSV to the directory the
manuscript reads. Scripts are idempotent: re-running them overwrites their own
outputs and nothing else.

---

## 5 Figure order

The authoritative figure order is `docs/manuscript_draft_v16.md`, which carries
the Markdown `![]()` embed next to each caption. The mapping from manuscript
figure number to file in `figures/` is:

| Figure | File | Figure | File |
|---|---|---|---|
| 1 | `fig0_roadmap.png` | 10 | `fig6_dependence_parameter.png` |
| 2 | `fig0_schematic.png` | 11 | `fig_rho_transfer.png` |
| 3 | `fig12_identified_set_geometry.png` | 12 | `fig12a_training_blocks.png` |
| 4 | `fig13_reference_windows.png` | 13 | `fig13_psd_asymptotic.png` |
| 5 | `fig1_scale_amplification.png` | 14 | `fig8_compositional.png` |
| 6 | `fig2_envelope_and_ladder.png` | 15 | `fig15_variogram_baseline.png` |
| 7 | `fig3_landscape_sweep.png` | 16 | `fig16_truncation_grading.png` |
| 8 | `fig4_calibration.png` | 17 | `figM1_range_ambiguity.png` |
| 9 | `fig5_real_landcover.png` | 18 | `figM2_anisotropy.png` |

---

## 6 Terminology

Two conventions are load-bearing and should be preserved by anyone reusing
this code:

1. **`normalised gap`, never `relative error`.** The quantity plotted in
   Figure 13b is the joint-mixability gap divided by the single-pixel marginal
   variance. It is an *absolute* gap normalised by a variance, not an error
   measured against the exact endpoint — which is undefined precisely in the
   regime where the relaxation returns zero while the sharp lower endpoint is
   positive.

2. **`enumerated reference value` for WorldCover, `truth` only for synthetic.**
   The synthetic landscapes have a constructed, error-free value: that is
   `truth`. WorldCover is a *classified reference field*, and the quantity
   enumerated from it is a spatial statistic of the mapped labels — an
   `enumerated reference value` or `reference-field value`. It is not ground
   truth, and it does not validate the posterior uncertainty of any
   probabilistic-composition model. The two must not share a name.

---

## 7 Licence and citation

Code released under the MIT Licence (see `LICENSE`). The manuscript text and
figures are copyright the authors.

If you use this code, please cite the manuscript (citation to be added on
publication) and ESA WorldCover 2021 v200 (Zanaga et al., 2022).
