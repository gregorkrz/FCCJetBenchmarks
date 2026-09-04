"""Trial 06 - four-param with a log(E) term: a/sqrt(E)+b+c/E+d*log(E).

Adds a slow logarithmic drift on top of the classic form to absorb residual
high-E curvature that the constant alone can't. 4 free params.

RESULTS: MAE=4.5401e-04 RMSE=5.6089e-04 maxAE=1.1059e-03 looCV_MAE=8.3772e-04 chi2/ndf=1.609e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    try:
        popt, _ = curve_fit(m_logE, x, y, p0=[0.5, 0.03, 0.1, 0.0], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("06_logE_4param",
                      "a/sqrt(E)+b+c/E+d*log(E), unweighted", fit_one)
