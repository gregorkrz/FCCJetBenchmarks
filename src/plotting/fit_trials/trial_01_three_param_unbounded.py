"""Trial 01 - three_param linear sum, UNBOUNDED, unweighted.

Baseline reference: a/sqrt(E) + b + c/E with no parameter bounds and ordinary
least squares. This is the honest "best RMSE" target the bounded production fit
should be compared against.

RESULTS: MAE=8.4170e-04 RMSE=1.0386e-03 maxAE=2.1871e-03 looCV_MAE=1.1685e-03 chi2/ndf=5.753e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("01_three_param_unbounded",
                      "a/sqrt(E)+b+c/E, unbounded, unweighted (LSQ)", fit_one)
