"""Trial 39 - a/sqrt(E)+b+c/E+e/E^1.5, weighted with 1% floor. 4 params.

The gentler-steep-term model (trial 29) with balanced weighting. A middle ground
between trial 20's aggressive 1/E^2 and the classic noise term, plus weighting
that keeps chi2 sane.

RESULTS: MAE=4.6382e-04 RMSE=6.2957e-04 maxAE=1.4264e-03 looCV_MAE=8.6151e-04 chi2/ndf=1.063e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_invE15_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_invE15_noise, x, y, p0=[0.5, 0.03, 0.1, 0.01],
                            sigma=sigma, absolute_sigma=True, maxfev=40000)
    except Exception:
        return None
    return lambda X: m_invE15_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("39_invE15_weighted_floor",
                      "a/sqrt(E)+b+c/E+e/E^1.5, weighted 1% floor", fit_one)
