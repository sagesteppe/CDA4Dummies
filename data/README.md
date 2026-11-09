# Data: 2026 ground-truth point predictions

`groundtruth_point_predictions.csv` — one row per (model x ground-truth point x AOA domain
x calibration). Source: `Ecoloradense/scripts/Second_groundtruth/evaluate_2026_pertree.R`
(per-point writeout), consolidated across three CV designs. This is the raw material for
decision-curve analysis by applicability domain and across assumed prevalence; it needs no
re-prediction (no models / HDD).

| column | meaning |
|---|---|
| `tree` | source results tree (raw design id) |
| `model` | model id: `<res>-Iteration<i>-PA1:<ratio>DO:0-Seed<seed>` |
| `design` | CV design: `classic` (results_classicsplit), `spatial` (spatialmapping split-tau), `lmpo` (LMPO split-tau) |
| `domain` | AOA domain: `Interp <=270m`, `Interp >270m`, `Extrapolation` |
| `point_id` | row index into the 2026 ground-truth point set (`gt`); join back for coordinates / patch keys |
| `Presence` | observed field truth, 1 = present, 0 = absent |
| `prob` | model predicted probability of occurrence |
| `dist_known_occ` | metres to nearest known occurrence (the <=270m / >270m split) |
| `DI` | AOA dissimilarity index of the point (continuous out-of-domain measure) |
| `DI_thresh` | model's AOA threshold; `DI > DI_thresh` => Extrapolation |
| `resolution` | predictor resolution (`1arc`, `3arc`, `1-3arc`) |
| `iteration` | modelling iteration |
| `SurfaceArea` / `Distance` | raw domain components `domain` is built from |
| `calibration` | `raw`, or `calibrated` (Platt, interpolation only) |

Notes:
- Extrapolation carries only `raw` (Platt is applied inside the AOA). Manuscript convention:
  calibrated interpolation + raw extrapolation.
- Prevalence here is the enriched field-sample prevalence (~0.34), NOT landscape prevalence.
  Reweight with case weights (presence = pi/n1, absence = (1-pi)/n0) to score net benefit at
  a target prevalence.
