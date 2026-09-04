"""Trial 30 - a/sqrt(E)+b+c/E^r  (free steep exponent r). 4 params.

Generalises trials 28/29: let the data pick the steep-term exponent r instead
of fixing it at 1, 1.5 or 2. Bounded r in [1,3] to keep it physically a
low-E-dominated term and avoid degeneracy with the sqrt term.

RESULTS: MAE=1.3817e-03 RMSE=1.7469e-03 maxAE=3.8181e-03 looCV_MAE=1.6222e-03 chi2/ndf=7.834e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_free_low_power


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    try:
        popt, _ = curve_fit(m_free_low_power, x, y, p0=[0.5, 0.03, 0.1, 1.5],
                            bounds=([0, 0, 0, 1.0], [np.inf, np.inf, np.inf, 3.0]), maxfev=40000)
    except Exception:
        return None
    return lambda X: m_free_low_power(X, *popt)


if __name__ == "__main__":
    harness.run_trial("30_free_low_power",
                      "a/sqrt(E)+b+c/E^r, free r in [1,3]", fit_one)
