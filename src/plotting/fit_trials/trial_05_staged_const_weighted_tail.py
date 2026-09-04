"""Trial 05 - STAGED with error-weighted tail estimate of the constant.

Like trial 03 but b is an inverse-variance-weighted mean of the last few points
(using the stored stat errors), which respects that the highest-E point is often
the lowest-N / noisiest. Then fix b and fit a,c weighted too.

RESULTS: MAE=3.5294e-03 RMSE=7.2149e-03 maxAE=2.1496e-02 looCV_MAE=1.9732e-03 chi2/ndf=1.507e+05 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness

N_TAIL = 4


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    e = err[order] if err is not None else None

    tail_y = y[-N_TAIL:]
    if e is not None and np.all(np.isfinite(e[-N_TAIL:])) and np.all(e[-N_TAIL:] > 0):
        w = 1.0 / e[-N_TAIL:] ** 2
        b = float(np.sum(w * tail_y) / np.sum(w))
    else:
        b = float(np.mean(tail_y))

    sigma = harness.safe_sigma(err[order] if err is not None else None, y, floor_frac=0.0)
    try:
        popt, _ = curve_fit(lambda E, a, c: a / np.sqrt(E) + b + c / E,
                            x, y, p0=[0.5, 0.1], sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    a, c = popt
    return lambda X: a / np.sqrt(np.asarray(X, float)) + b + c / np.asarray(X, float)


if __name__ == "__main__":
    harness.run_trial("05_staged_const_weighted_tail",
                      f"inverse-variance b from last {N_TAIL} pts, then weighted a,c fit", fit_one)
