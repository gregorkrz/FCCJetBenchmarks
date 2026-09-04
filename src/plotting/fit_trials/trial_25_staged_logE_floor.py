"""Trial 25 - STAGED constant + log(E) model, weighted with 1% floor.

The kitchen-sink of the promising ideas:
  1. Pin b from the inverse-variance-weighted tail (trial 05 idea).
  2. Seed a,c,d for a/sqrt(E)+b+c/E+d*log(E) with b fixed.
  3. Release all 4 and refit, weighted with a 1% error floor.

RESULTS: MAE=5.2928e-04 RMSE=7.4436e-04 maxAE=1.8162e-03 looCV_MAE=8.4751e-04 chi2/ndf=1.214e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE

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
        b0 = float(np.sum(w * tail_y) / np.sum(w))
    else:
        b0 = float(np.mean(tail_y))

    sigma = harness.safe_sigma(e, y, floor_frac=0.01)
    try:
        # seed with b fixed
        p_seed, _ = curve_fit(lambda E, a, c, d: a / np.sqrt(E) + b0 + c / E + d * np.log(E),
                            x, y, p0=[0.5, 0.1, 0.0], sigma=sigma, absolute_sigma=True, maxfev=20000)
        a0, c0, d0 = p_seed
        popt, _ = curve_fit(m_logE, x, y, p0=[a0, b0, c0, d0],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("25_staged_logE_floor",
                      "staged b + log(E) 4-param, weighted, 1% floor", fit_one)
