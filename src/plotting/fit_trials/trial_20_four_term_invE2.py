"""Trial 20 - a/sqrt(E)+b+c/E+e/E^2 (extra steep 1/E^2 low-E term).

Adds an even steeper term than noise to better capture the sharp low-E turn-up
that the classic form tends to undershoot. 4 params.

RESULTS: MAE=3.9072e-04 RMSE=4.9512e-04 maxAE=1.0172e-03 looCV_MAE=8.2180e-04 chi2/ndf=1.176e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    try:
        popt, _ = curve_fit(m_stoch_const_noise2, x, y, p0=[0.5, 0.03, 0.1, 0.01], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("20_four_term_invE2",
                      "a/sqrt(E)+b+c/E+e/E^2", fit_one)
