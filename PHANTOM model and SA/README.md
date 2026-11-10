# Tide-mediated propagule dispersal and establishment model: reproduction package

From the raw Gazi Bay survey workbook to the figures of the paper, in one command.

## Requirements

- Python 3.10+ with numpy, pandas, scipy, matplotlib, pillow, xlrd
- R 4.x, only if the adult stand is rebuilt (set the environment variable `RSCRIPT` to your `Rscript` executable)

## Run

```bash
python run_all.py                     # all model runs and figures
python run_all.py --rebuild-inputs    # first rebuild the shipped inputs from the raw data (needs R, about 15 min)
python run_all.py --figures-only      # redraw the figures from existing runs
```

`run_all.py` makes 128 one-year model runs (about 5 minutes per run per worker, `--workers 12` by default)
and writes about 6 GB of run tables to `results/`. The figures appear in `results/figures/`.

A single model run:

```bash
python tidal_dispersal_model.py -s a --days 365 --trees results/reconstruction/stand3T_T4T3T1T6_single/demo_reconstructed_stand3T.csv
```

## Contents, by section of the methods

| methods | what | files |
|---|---|---|
| raw data | survey workbook: adult census, 1997 and 1999 juvenile surveys | `GAZ1_6 (1).XLS` |
| raw data | published map (Dahdouh-Guebas et al. 2004, Fig. 6); hand-drawn forest edge | `results/reconstruction/figures/dahdouh_guebas_2004_page09_300dpi.png`, `reconstruction/forest_edge_user.csv` |
| juvenile data | workbook -> quadrat tables `juvenile_*_t1_5.csv`, zoned by position | `extract_juvenile.py`, `zone_quadrats_by_position.py` |
| adult stand | PCQM records of transects T1, T3, T4, T6 (T2 withheld) -> reconstructed stand, five members | `reconstruction/run_stand3T.py` and the modules beside it; substrate map from `reconstruction/make_substrate_map.py` |
| 2.3 model and parameters | the model with every parameter constant | `tidal_dispersal_model.py`, `docs/parameter_provenance.md` |
| 2.4.1 drag sensitivity | five drag profiles x 3 seeds | `drag_sensitivity.py` -> `figure_drag_zonation.py` |
| 2.4.2 deconstruction | species-specific vs equalised transport traits, three origin groups | `transport_treatments.py` -> `figure_endpoints_main.py` |
| 2.4.2 one-at-a-time sweep | draft, mobility, lateral spread, ten contrast levels | `equal_cv_experiment.py` (runs and figure) |
| 2.4.2 Sonneratia curves | supply x 1000, distribution curves only | `run_sonneratia_shape.py` |
| establishment | stranded vs established per species | `figure_establishment_filter.py` |
| evaluation | 1999 and 1997 surveys vs model; withheld adults vs recruits | `virtual_survey.py`, `evaluate_recruits_vs_adults.py` -> `figure_evaluation_stacked.py`, `figure_evaluation_1999_change_adults.py` |
| study site | published map next to the reconstructed stand | `fig_layout_stand.py` |

Supporting modules imported by the scripts above: `transect_geometry.py` (transect layout, survey
operator), `cascade_metrics.py` and `validate_1999_patterns.py` (per-run summaries). In `reconstruction/`:
`prep_pcqm.py`, `run_mtppr.py`, `mtppr_tile.R`, `postprocess.py`, `recon_config.py`, `recon_geometry.py`,
`virtual_pcqm.py`. The point-pattern engine in `ChrisWudel-Multi-trait-point-pattern-reconstruction-6982c39/`
is the unmodified MTPPR v2.0.0 of Chris Wudel (licence included); its checksum is asserted before every run.

Shipped derived inputs, all rebuilt by `--rebuild-inputs`: the juvenile tables, `results/reconstruction/substrate_map.npz`
and the stand in `results/reconstruction/stand3T_T4T3T1T6_single/` (member 2, `demo_reconstructed_stand3T.csv`,
2,119 trees, is the stand the model runs on).

## Reproducibility

- Model runs are deterministic for a given seed. Seeds: 42, 49, 56 for the experiments; 20 seeds
  (42 to 175 in steps of 7) for the calibrated model.
- `--rebuild-inputs` reproduces the shipped juvenile tables, substrate map and all five stand members
  byte for byte (checked on Windows, R 4.5.3).
- The Sonneratia runs with supply x 1000 are used only to draw its distribution curves, never for
  densities or recruit fractions.
- The calibration runs that selected the mobility values, the boundary supply and the Sonneratia supply
  multiplier are not part of this package; the calibrated values are the constants in `tidal_dispersal_model.py`.
