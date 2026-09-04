"""Trial 28 - a/sqrt(E)+b+c/E^1.5  (single steeper term, still 3 params).

Round 1 showed the low-E turn-up needs a steeper term than 1/E. Instead of
adding a 4th param, just make the steep term 1/E^1.5. Same param count as the
production three_param but should track the turn-up better.

RESULTS: MAE=8.0750e-04 RMSE=9.9237e-04 maxAE=2.0280e-03 looCV_MAE=1.1607e-03 chi2/ndf=5.269e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_invE15


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    try:
        popt, _ = curve_fit(m_invE15, x, y, p0=[0.5, 0.03, 0.1], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_invE15(X, *popt)


if __name__ == "__main__":
    harness.run_trial("28_invE15",
                      "a/sqrt(E)+b+c/E^1.5, 3 params", fit_one)
