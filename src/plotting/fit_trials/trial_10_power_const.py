"""Trial 10 - a*E^(-p) + b  (free exponent instead of fixed 1/sqrt(E)).

Lets the data choose the energy-scaling exponent p (calorimeter stochastic term
is nominally p=0.5, but jet-level effects can shift it). 3 params.

RESULTS: MAE=1.3504e-03 RMSE=1.6902e-03 maxAE=3.6787e-03 looCV_MAE=1.5945e-03 chi2/ndf=8.109e+04 (fit 159, failed 0)
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
    try:
        popt, _ = curve_fit(m_power_const, x, y, p0=[0.5, 0.5, 0.03],
                            bounds=([0, 0.1, 0], [np.inf, 2.0, np.inf]), maxfev=20000)
    except Exception:
        return None
    return lambda X: m_power_const(X, *popt)


if __name__ == "__main__":
    harness.run_trial("10_power_const",
                      "a*E^(-p)+b, free exponent p in [0.1,2]", fit_one)
