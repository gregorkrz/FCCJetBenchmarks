"""Trial 26 - a/sqrt(E)+b+c/E+e/E^2+d*log(E)  (BOTH round-1 winners, 5 params).

Round 1's two best ideas were the steep e/E^2 term (trial 20, best MAE) and the
d*log(E) drift (trial 06). Combine them. 5 params - safe for the >=9-point
series (ndf>=4) but watch for overfit on the shortest ones.

RESULTS: MAE=2.5478e-04 RMSE=3.2658e-04 maxAE=6.7990e-04 looCV_MAE=8.9428e-04 chi2/ndf=6.454e+03 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE_and_invE2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 6:
        return None
    try:
        popt, _ = curve_fit(m_logE_and_invE2, x, y, p0=[0.5, 0.03, 0.1, 0.01, 0.0], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_logE_and_invE2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("26_logE_and_invE2",
                      "a/sqrt(E)+b+c/E+e/E^2+d*log(E), 5 params", fit_one)
