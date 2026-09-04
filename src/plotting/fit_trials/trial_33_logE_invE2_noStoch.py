"""Trial 33 - b + c/E + e/E^2 + d*log(E)  (drop the 1/sqrt(E) term). 4 params.

Round-1 fits often drove the stochastic a-term toward small values. Test whether
the sqrt term is redundant once steep (c/E, e/E^2) and slow (log) terms are
present - i.e. is a/sqrt(E) actually earning its parameter?

RESULTS: MAE=4.0570e-04 RMSE=5.1114e-04 maxAE=1.0371e-03 looCV_MAE=8.2960e-04 chi2/ndf=1.243e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE_invE2_noStoch


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    try:
        popt, _ = curve_fit(m_logE_invE2_noStoch, x, y, p0=[0.03, 0.5, 0.05, 0.0], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_logE_invE2_noStoch(X, *popt)


if __name__ == "__main__":
    harness.run_trial("33_logE_invE2_noStoch",
                      "b+c/E+e/E^2+d*log(E), no sqrt term", fit_one)
