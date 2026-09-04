"""Trial 13 - quadrature sum sqrt((a/sqrt(E))^2 + b^2 + (c/E)^2).

The physically standard way calorimeter resolution terms combine (added in
quadrature, not linearly). Same 3 params as the classic linear form.

RESULTS: MAE=1.6397e-03 RMSE=2.0977e-03 maxAE=4.5529e-03 looCV_MAE=1.8200e-03 chi2/ndf=1.396e+05 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_quad3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    try:
        popt, _ = curve_fit(m_quad3, x, y, p0=[0.5, 0.03, 0.1], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_quad3(X, *popt)


if __name__ == "__main__":
    harness.run_trial("13_quad3",
                      "sqrt((a/sqrt(E))^2+b^2+(c/E)^2), quadrature, unweighted", fit_one)
