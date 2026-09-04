"""Trial 27 - trial-20 model (e/E^2), weighted with 1% error floor.

Round 1's outright winner (trial 20, a/sqrt(E)+b+c/E+e/E^2) was unweighted.
Add balanced weighting (1% floor) to see if it improves the statistically
meaningful chi2 while keeping the excellent point-wise error.

RESULTS: MAE=4.4943e-04 RMSE=6.0452e-04 maxAE=1.3485e-03 looCV_MAE=8.7698e-04 chi2/ndf=1.049e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_stoch_const_noise2, x, y, p0=[0.5, 0.03, 0.1, 0.01],
                            sigma=sigma, absolute_sigma=True, maxfev=40000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("27_invE2_weighted_floor",
                      "a/sqrt(E)+b+c/E+e/E^2, weighted 1% floor", fit_one)
