"""Trial 21 - STAGED 3-way: pin b from tail, pin c from low-E, then fit a.

A more aggressive staged decomposition:
  1. b = mean of highest-E points (plateau).
  2. Subtract b; on the lowest-E points, sigma/E - b ~ c/E (noise dominates),
     so estimate c from (y-b)*E at low E.
  3. Fix b and c, fit only the stochastic a on all points.
Then one final light release of all three from that seed.

RESULTS: MAE=8.4170e-04 RMSE=1.0386e-03 maxAE=2.1871e-03 looCV_MAE=1.1685e-03 chi2/ndf=5.753e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise

N_TAIL = 3
N_HEAD = 3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    b0 = float(np.mean(y[-N_TAIL:]))
    c0 = float(np.mean((y[:N_HEAD] - b0) * x[:N_HEAD]))
    c0 = max(c0, 0.0)
    try:
        a_fit, _ = curve_fit(lambda E, a: a / np.sqrt(E) + b0 + c0 / E,
                            x, y, p0=[0.5], maxfev=20000)
        a0 = float(a_fit[0])
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[a0, b0, c0], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("21_staged_const_noise_then_stoch",
                      "pin b (tail) & c (head), fit a, then release all", fit_one)
