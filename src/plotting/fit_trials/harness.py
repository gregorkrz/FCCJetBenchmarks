"""Shared harness for JER fit experiments (see fit_trials/README or the goal).

Each trial script (trial_XX_*.py) imports this, defines a `fit_one(x, y, err)`
that returns a prediction callable `pred(x)` (or None on failure), and calls
`run_trial(NAME, DESC, fit_one)`. The harness handles: loading every energy-JER
(mid_points, sigma_over_E, per-point error) series from dashboard_data.json,
evaluating the trial on all of them, aggregating error metrics, printing a
report, and appending one row to fit_trials/results.jsonl so all trials can be
ranked together afterwards.

Design goals:
  * Trial scripts stay tiny and self-documenting - just the fitting idea.
  * Metrics are computed identically across trials (fair comparison).
  * Runnable standalone:  python src/plotting/fit_trials/trial_XX_name.py
"""
import json
import os
import sys

import numpy as np

# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

_DEFAULT_JSON = os.environ.get(
    "DASHBOARD_JSON",
    "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks_histograms/IDEA_20260120/plots/dashboard_data.json",
)

RESULTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.jsonl")

_CACHE = {}


def load_series(json_path=None):
    """Return a list of dicts: {method, process, x, y, err} for every energy._all.

    x = mid_points (GeV), y = sigma/E, err = sigma_over_E_err (1-sigma per point,
    may contain zeros/NaN where a bin had no error). Cached per path.
    """
    json_path = json_path or _DEFAULT_JSON
    if json_path in _CACHE:
        return _CACHE[json_path]
    with open(json_path) as fh:
        d = json.load(fh)
    series = []
    for method, mv in d["methods"].items():
        for process, pv in mv["processes"].items():
            entry = (pv.get("energy") or {}).get("_all")
            if not entry:
                continue
            x = np.asarray(entry.get("mid_points", []), dtype=float)
            y = np.asarray(entry.get("sigma_over_E", []), dtype=float)
            if x.size < 3:
                continue
            err = entry.get("sigma_over_E_err")
            err = np.asarray(err, dtype=float) if err is not None else None
            series.append({"method": method, "process": process, "x": x, "y": y, "err": err})
    _CACHE[json_path] = series
    return series


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------

def metrics_for(x, y, err, pred):
    """Error metrics for one series given a prediction callable pred(x)->y_hat."""
    yhat = np.asarray(pred(x), dtype=float)
    resid = y - yhat
    n = resid.size
    out = {
        "MAE": float(np.mean(np.abs(resid))),
        "MSE": float(np.mean(resid ** 2)),
        "RMSE": float(np.sqrt(np.mean(resid ** 2))),
        "maxAE": float(np.max(np.abs(resid))),
        "n": int(n),
    }
    if err is not None and err.shape == resid.shape:
        good = np.isfinite(err) & (err > 0)
        ndf = int(good.sum()) - 3  # nominal 3 free params; comparative only
        if ndf > 0:
            out["chi2_ndf"] = float(np.sum((resid[good] / err[good]) ** 2) / ndf)
    return out


def _aggregate(rows):
    """Mean of each metric across all series (skipping failed fits)."""
    keys = ["MAE", "MSE", "RMSE", "maxAE"]
    agg = {}
    for k in keys:
        vals = [r[k] for r in rows if r is not None]
        agg[k] = float(np.mean(vals)) if vals else None
    chi = [r["chi2_ndf"] for r in rows if r is not None and "chi2_ndf" in r]
    agg["chi2_ndf"] = float(np.mean(chi)) if chi else None
    agg["n_series"] = sum(1 for r in rows if r is not None)
    agg["n_failed"] = sum(1 for r in rows if r is None)
    return agg


# --------------------------------------------------------------------------
# Trial runner
# --------------------------------------------------------------------------

def loo_cv_mae(x, y, err, fit_one):
    """Leave-one-out CV MAE: refit on all-but-one point, predict the held-out one.

    Guards against the in-sample MAE being flattered by high parameter counts on
    short (9-point) series. Returns (loo_mae, n_folds) or (None, 0) if too short
    or fits fail. Interior points only are held out (endpoints extrapolate wildly
    for these steep models and would unfairly punish every model equally).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    if n < 6:
        return None, 0
    order = np.argsort(x)
    x, y = x[order], y[order]
    e = np.asarray(err, float)[order] if err is not None else None
    abs_errs = []
    for i in range(1, n - 1):  # skip the two endpoints
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        ek = e[keep] if e is not None else None
        try:
            pred = fit_one(x[keep], y[keep], ek)
        except Exception:
            pred = None
        if pred is None:
            continue
        try:
            yhat = float(np.asarray(pred(np.array([x[i]])))[0])
        except Exception:
            continue
        if np.isfinite(yhat):
            abs_errs.append(abs(y[i] - yhat))
    if not abs_errs:
        return None, 0
    return float(np.mean(abs_errs)), len(abs_errs)


def run_trial(name, desc, fit_one, json_path=None, verbose=True, append=True, do_loo=True):
    """Evaluate `fit_one` over every series and report + persist aggregate metrics.

    fit_one(x, y, err) -> pred callable (pred(x)->yhat) or None on failure.
    Returns the aggregate metrics dict. When do_loo is set, also computes a
    leave-one-out CV MAE (mean over series) as an overfitting guard.
    """
    series = load_series(json_path)
    per = []
    loo_vals = []
    for s in series:
        try:
            pred = fit_one(s["x"], s["y"], s["err"])
        except Exception:
            pred = None
        if pred is None:
            per.append(None)
        else:
            try:
                per.append(metrics_for(s["x"], s["y"], s["err"], pred))
            except Exception:
                per.append(None)
        if do_loo:
            lmae, nf = loo_cv_mae(s["x"], s["y"], s["err"], fit_one)
            if lmae is not None:
                loo_vals.append(lmae)
    agg = _aggregate(per)
    agg["trial"] = name
    agg["desc"] = desc
    agg["loo_mae"] = float(np.mean(loo_vals)) if loo_vals else None

    if verbose:
        print(f"=== {name} ===")
        print(desc)
        print(f"  series fit: {agg['n_series']}   failed: {agg['n_failed']}")
        print(f"  MAE   = {agg['MAE']:.4e}" if agg["MAE"] is not None else "  MAE   = n/a")
        print(f"  RMSE  = {agg['RMSE']:.4e}" if agg["RMSE"] is not None else "  RMSE  = n/a")
        print(f"  MSE   = {agg['MSE']:.4e}" if agg["MSE"] is not None else "  MSE   = n/a")
        print(f"  maxAE = {agg['maxAE']:.4e}" if agg["maxAE"] is not None else "  maxAE = n/a")
        chi = agg.get("chi2_ndf")
        print(f"  chi2/ndf(mean) = {chi:.4e}" if chi is not None else "  chi2/ndf = n/a")

    if append:
        # De-dupe: drop any prior row for this trial name, then append fresh.
        rows = []
        if os.path.exists(RESULTS_PATH):
            with open(RESULTS_PATH) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        r = json.loads(line)
                    except Exception:
                        continue
                    if r.get("trial") != name:
                        rows.append(r)
        rows.append(agg)
        with open(RESULTS_PATH, "w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")

    return agg


# --------------------------------------------------------------------------
# Shared fitting helpers (used by many trials)
# --------------------------------------------------------------------------

def clean(x, y, err):
    """Drop non-finite points; return (x, y, err_or_None) with err made safe.

    err entries that are non-finite or <=0 are replaced by NaN so callers can
    decide how to handle them (curve_fit rejects sigma<=0).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if err is None:
        return x, y, None
    err = np.asarray(err, float)[m]
    err = np.where(np.isfinite(err) & (err > 0), err, np.nan)
    return x, y, err


def safe_sigma(err, y, floor_frac=0.0):
    """Turn a raw err array into a curve_fit sigma with an optional relative floor.

    sigma_i = sqrt(err_i^2 + (floor_frac * y_i)^2); NaN errs fall back to the
    floor (or to 1.0 if no floor). Returns None if all sigmas are unusable.
    """
    if err is None and floor_frac <= 0:
        return None
    y = np.asarray(y, float)
    base = np.where(np.isfinite(err), err, 0.0) if err is not None else np.zeros_like(y)
    floor = floor_frac * np.abs(y)
    sigma = np.sqrt(base ** 2 + floor ** 2)
    sigma = np.where(sigma > 0, sigma, np.nan)
    if not np.any(np.isfinite(sigma)):
        return None
    # Replace remaining NaN with the max finite sigma (least informative).
    fill = np.nanmax(sigma) if np.any(np.isfinite(sigma)) else 1.0
    return np.where(np.isfinite(sigma), sigma, fill)
