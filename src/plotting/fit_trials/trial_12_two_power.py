"""Trial 12 - a*E^(-p) + c*E^(-q)  (two free power laws, no constant).

A purely empirical two-term power law. No constant floor - tests whether these
data even need a non-vanishing high-E plateau, or if two decaying terms suffice.

RESULTS: MAE=1.2928e-03 RMSE=1.6214e-03 maxAE=3.5494e-03 looCV_MAE=1.5465e-03 chi2/ndf=7.223e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_two_power


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    try:
        popt, _ = curve_fit(m_two_power, x, y, p0=[0.5, 0.5, 0.1, 1.0],
                            bounds=([0, 0.05, 0, 0.05], [np.inf, 2.0, np.inf, 3.0]), maxfev=30000)
    except Exception:
        return None
    return lambda X: m_two_power(X, *popt)


if __name__ == "__main__":
    harness.run_trial("12_two_power",
                      "a*E^(-p)+c*E^(-q), two free power laws", fit_one)
