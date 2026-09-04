"""Trial 07 - a/sqrt(E)+b+d*log(E)  (log replaces the 1/E noise term).

Tests whether a log(E) drift is a better 3rd degree of freedom than the steep
c/E noise term - same param count as the classic three_param.

RESULTS: MAE=9.4390e-04 RMSE=1.1830e-03 maxAE=2.6751e-03 looCV_MAE=1.2307e-03 chi2/ndf=6.923e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_logE


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    try:
        popt, _ = curve_fit(m_stoch_const_logE, x, y, p0=[0.5, 0.03, 0.0], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("07_logE_noNoise",
                      "a/sqrt(E)+b+d*log(E) (log instead of c/E)", fit_one)
