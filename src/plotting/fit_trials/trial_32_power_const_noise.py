"""Trial 32 - a*E^(-p) + b + c/E  (free stochastic exponent + constant + noise).

Round 1's free-exponent power model (trial 10) lacked a noise term. Add c/E back
so the low-E turn-up is captured by an explicit steep term while the main
exponent p stays free. 4 params.

RESULTS: MAE=1.3493e-03 RMSE=1.6891e-03 maxAE=3.6781e-03 looCV_MAE=1.5937e-03 chi2/ndf=8.104e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_power_const_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    try:
        popt, _ = curve_fit(m_power_const_noise, x, y, p0=[0.5, 0.5, 0.03, 0.1],
                            bounds=([0, 0.1, 0, 0], [np.inf, 2.0, np.inf, np.inf]), maxfev=40000)
    except Exception:
        return None
    return lambda X: m_power_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("32_power_const_noise",
                      "a*E^(-p)+b+c/E, free stochastic exponent", fit_one)
