"""Trial 08 - a/sqrt(E)+b+d*log(E)/E  (log-suppressed noise term).

A physically softer noise term: log(E)/E decays a touch slower than 1/E at low E
but still vanishes at high E, unlike a bare log which grows.

RESULTS: MAE=8.7399e-04 RMSE=1.0669e-03 maxAE=2.1934e-03 looCV_MAE=1.3645e-03 chi2/ndf=5.071e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_inv_sqrt_log


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    try:
        popt, _ = curve_fit(m_inv_sqrt_log, x, y, p0=[0.5, 0.03, 0.1], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_inv_sqrt_log(X, *popt)


if __name__ == "__main__":
    harness.run_trial("08_logE_over_E",
                      "a/sqrt(E)+b+d*log(E)/E", fit_one)
