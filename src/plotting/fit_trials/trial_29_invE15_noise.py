"""Trial 29 - a/sqrt(E)+b+c/E+e/E^1.5  (noise + a milder steep term). 4 params.

Between trial 20 (e/E^2, very steep) and the classic form: a 1/E^1.5 term is a
gentler extra low-E degree of freedom that may generalise better than 1/E^2.

RESULTS: MAE=3.9955e-04 RMSE=5.0420e-04 maxAE=1.0279e-03 looCV_MAE=8.1079e-04 chi2/ndf=1.231e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_invE15_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    try:
        popt, _ = curve_fit(m_invE15_noise, x, y, p0=[0.5, 0.03, 0.1, 0.01], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_invE15_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("29_invE15_noise",
                      "a/sqrt(E)+b+c/E+e/E^1.5, 4 params", fit_one)
