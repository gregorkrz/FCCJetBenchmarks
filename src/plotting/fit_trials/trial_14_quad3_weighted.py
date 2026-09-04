"""Trial 14 - quadrature 3-term, error-weighted.

Trial 13's quadrature model with stat-error weighting.

RESULTS: MAE=1.9176e-03 RMSE=4.1817e-03 maxAE=1.2575e-02 looCV_MAE=1.0673e-03 chi2/ndf=3.452e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_quad3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.0)
    try:
        popt, _ = curve_fit(m_quad3, x, y, p0=[0.5, 0.03, 0.1],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_quad3(X, *popt)


if __name__ == "__main__":
    harness.run_trial("14_quad3_weighted",
                      "quadrature 3-term, stat-error weighted", fit_one)
