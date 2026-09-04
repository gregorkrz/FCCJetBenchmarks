# JER fit experiments

Experiments to find a better energy-dependence fit for the jet energy resolution
(sigma_E/E vs E) than the production `three_param` `A/sqrt(E) + B + C/E` with
`B` bounded to `[0.005, 0.04]`.

## Layout

- `harness.py` - loads every energy-JER (mid_points, sigma/E, per-point error)
  series from `dashboard_data.json`, and provides `run_trial()`, the error
  metrics (MAE/MSE/RMSE/maxAE, chi2/ndf), and a **leave-one-out CV MAE**
  (`loo_mae`) that refits on all-but-one interior point to catch overfitting.
- `models.py` - the candidate model functions (physics letters, not the swapped
  A/B/C display labels).
- `trial_XX_*.py` - one strategy each; runnable standalone
  (`python trial_XX_*.py`) and self-documenting. Its headline metrics are
  back-filled into the `RESULTS:` line of its own docstring.
- `run_all.py` - runs every trial, writes `results.jsonl`, prints a leaderboard.
  `--metric loo_mae` ranks by the overfitting-robust CV metric (recommended).

## Reproduce

```bash
source env.sh   # sets DASHBOARD_JSON target via PATH_TO_HISTOGRAMS
cd src/plotting/fit_trials
python3 run_all.py --metric loo_mae     # cross-validated ranking
python3 run_all.py --metric MAE         # in-sample ranking (can mislead!)
```

Point a run at a different dataset with `DASHBOARD_JSON=/path/to/dashboard_data.json`.

## Headline result (159 method x process series)

Ranked by **leave-one-out CV MAE** (generalization, not memorization):

| Rank | Trial | Model | CV MAE | in-sample MAE |
|------|-------|-------|--------|---------------|
| 1 | `09_logE_4param_weighted` | `A/sqrt(E)+B+C/E+D*log(E)`, error-weighted | **6.5e-4** | 6.8e-4 |
| 2 | `29_invE15_noise` | `A/sqrt(E)+B+C/E+e/E^1.5` | 8.1e-4 | 4.0e-4 |
| 3 | `20_four_term_invE2` | `A/sqrt(E)+B+C/E+e/E^2` | 8.2e-4 | 3.9e-4 |
| ... | | | | |
| 43 | `45_three_param_baseline_bounded` | **production incumbent** | 2.0e-3 | 1.8e-3 |

### Lessons

1. **A `d*log(E)` term is the best single added degree of freedom** and, when
   error-weighted (trial 09), generalizes best - CV error ~= in-sample error,
   3.1x better than the production fit on held-out points.
2. **In-sample MAE is misleading.** The 5-param `log + e/E^2` model (trial 26)
   had the best in-sample MAE (2.5e-4) but overfits: its CV MAE (8.9e-4) is
   nearly 3x worse. Always rank by `loo_mae`.
3. **The production bounds hurt.** Unbounded `three_param` (trial 01) alone
   already beats the bounded incumbent (trial 45) ~2x.
4. **Staged constant-first (the fix-b-from-the-tail idea):** fixing `b` from the
   noisy tail and never releasing it hurts (trials 03/22); using it only as a
   *seed* then releasing all params recovers the unbounded optimum (trial 04).
