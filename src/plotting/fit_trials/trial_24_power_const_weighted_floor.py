"""Trial 24 - a*E^(-p)+b, weighted with 1% error floor.

Combines the free-exponent power model (trial 10) with balanced weighting
(trial 15's 1% floor). Tests whether a data-driven exponent plus sane weighting
is the sweet spot.

RESULTS: MAE=1.3632e-03 RMSE=2.1271e-03 maxAE=5.2879e-03 looCV_MAE=1.3012e-03 chi2/ndf=4.908e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_power_const


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_power_const, x, y, p0=[0.5, 0.5, 0.03],
                            bounds=([0, 0.1, 0], [np.inf, 2.0, np.inf]),
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_power_const(X, *popt)


if __name__ == "__main__":
    harness.run_trial("24_power_const_weighted_floor",
                      "a*E^(-p)+b, weighted, 1% error floor", fit_one)
