"""Trial 31 - trial-20 model with a staged constant seed, then full release.

Marries the round-1 winner (e/E^2 term) with the user's staged-constant idea
done RIGHT (as a seed only, per trial 04's lesson): pin b from the tail to seed,
then release all 4 params. Tests whether better initialisation helps the 4-param
fit converge to a lower minimum.

RESULTS: MAE=3.9072e-04 RMSE=4.9512e-04 maxAE=1.0172e-03 looCV_MAE=8.2180e-04 chi2/ndf=1.176e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2

N_TAIL = 3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    b0 = float(np.mean(ys[-N_TAIL:]))
    try:
        # seed a,c,e with b fixed
        p_seed, _ = curve_fit(
            lambda E, a, c, e: a / np.sqrt(E) + b0 + c / E + e / E ** 2,
            xs, ys, p0=[0.5, 0.1, 0.01], maxfev=40000)
        a0, c0, e0 = p_seed
        popt, _ = curve_fit(m_stoch_const_noise2, xs, ys, p0=[a0, b0, c0, e0], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("31_invE2_seed_release",
                      "e/E^2 model, staged-b seed then release all 4", fit_one)
