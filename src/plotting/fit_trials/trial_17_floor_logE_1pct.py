"""Trial 17 - 4-param log(E) model, weighted with 1% error floor.

Combines the two most promising ingredients so far: the extra log(E) DOF (trial
06) and the balanced weighting of a 1% error floor (trial 15).

RESULTS: MAE=5.2928e-04 RMSE=7.4436e-04 maxAE=1.8162e-03 looCV_MAE=8.4751e-04 chi2/ndf=1.214e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_logE, x, y, p0=[0.5, 0.03, 0.1, 0.0],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("17_floor_logE_1pct",
                      "a/sqrt(E)+b+c/E+d*log(E), weighted, 1% error floor", fit_one)
