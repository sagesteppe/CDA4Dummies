# Example model — one SDM, for the "a priori vs realized" DCA example

A single worked example SDM for *Eriogonum coloradense*, pulled in so the Narrative can show
the first (pre-ground-truth) DCA formulation and compare it to what ground-truthing revealed.

- **Model id:** `1arc-Iteration1-PA1:1.023DO:0-Seed5101`
- **Source:** `Ecoloradense/Results_LMPO_splittau_phase2/tables/` (files renamed here for clarity;
  originals keep the model id in the filename).
- **PA ratio:** `PA1:1.023` = 1 presence : 1.023 absences in the training/test design, so the
  **test-holdout prevalence ≈ 1 / (1 + 1.023) = 0.4943** (near-balanced).
- **Landscape prevalence for this example:** **π = 0.025** (set by RCB for this portion).

| file | one row per | key fields |
|---|---|---|
| `test_holdout_metrics.csv` | metric | `pr_auc` (0.836), `roc_auc` (0.789), `brier_class`, `boyce`, `cox_intercept/slope`. Measured at the **test** prevalence (~0.494). This is the number reported *before* ground-truthing. |
| `groundtruth_metrics_by_domain.csv` | metric × AOA domain | realized `pr_auc`, `roc_auc`, `prevalence`, cox calibration, and pre-computed `dca_*` (prop_useful, mean_excess_nb, nb at pt .25/.50/.75). Domains: `Interpolation <=270m` (n=58), `Interpolation >270m` (n=178), `Extrapolation` (n=47). |
| `groundtruth_points.csv` | ground-truth point | `Presence`, `prob`, `dist_known_occ`, `DI`, `DI_thresh`, domain. Point-level material for this one model. |

## Why test PR-AUC must be corrected to landscape prevalence

Recall (TPR) and FPR are prevalence-invariant; **precision is not**. The test PR-AUC (0.836) is
read at ~0.494 prevalence, and at a landscape prevalence of 0.025 Whereby the model's precision — and
therefore its net benefit, is much lower. Realized in-domain PR-AUC here is 0.60, i.e. the
reported number is optimistic even before the prevalence correction. Correcting precision to a
target prevalence π:

    precision_π = (π · recall) / (π · recall + (1 − π) · FPR)

An exact re-derivation of the whole PR curve at a new π needs the test-holdout ROC/PR curve
(the per-point test predictions), not just the scalar AP — those live in the modeling step, not
in these files. The scalars here support the point-wise correction and the headline comparison.
