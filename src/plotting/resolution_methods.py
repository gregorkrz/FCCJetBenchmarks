"""Pluggable resolution-extraction and energy-dependence fit models.

This module is the single place where "how do we turn a histogram into a
resolution number" (SIGMA_METHODS) and "how does the resolution depend on
energy/angle" (RESOLUTION_MODELS) are defined. Both are plain dict
registries so that new algorithms can be added without touching the
plotting code in resolution_plots.py.

To add a new per-bin resolution estimator:
    def sigma_my_method(y, edges, wmin=0.7, wmax=1.2, percentage=0.683, epsilon=0.0005):
        ...
        return sigma, low, high, mpv   # or None if extraction failed
    SIGMA_METHODS["my_method"] = sigma_my_method

To add a new energy-dependence functional form:
    def my_model(E, a, b):
        return a * E + b
    RESOLUTION_MODELS["my_model"] = dict(func=my_model, p0=[1.0, 0.0],
                                          bounds=([-np.inf, -np.inf], [np.inf, np.inf]))

Then pass sigma_method="my_method" / model="my_model" to the plotting code.
"""
import numpy as np
from scipy.optimize import curve_fit

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def decorator(func):
            return func

        return decorator


@njit(cache=True)
def _find_narrowest_interval(centers, cumulative_weights, percentage, epsilon, wmin, wmax):
    """
    Find the narrowest interval containing the specified percentage of the distribution.
    Returns: (width, low, high) or (100.0, wmin, wmax) if no valid interval found.
    """
    n = len(centers)
    best_width = 100.0
    best_low = wmin
    best_high = wmax

    for i in range(n):
        for j in range(i, n):
            wy = cumulative_weights[j] - cumulative_weights[i]
            if abs(wy - percentage) < epsilon:
                wx = centers[j] - centers[i]
                if wx < best_width:
                    best_width = wx
                    best_low = centers[i]
                    best_high = centers[j]

    return best_width, best_low, best_high


def _narrowest_interval_fallback(centers, cumulative, percentage, epsilon, wmin, wmax):
    width = 100.0
    low, high = wmin, wmax
    for i in range(len(centers)):
        for j in range(i, len(centers)):
            wy = cumulative[j] - cumulative[i]
            if abs(wy - percentage) < epsilon:
                wx = centers[j] - centers[i]
                if wx < width:
                    low, high, width = centers[i], centers[j], wx
    return width, low, high


def sigma_std68(y, edges, wmin=0.7, wmax=1.2, percentage=0.683, epsilon=None):
    """Narrowest interval containing `percentage` of the distribution."""
    if epsilon is None:
        epsilon = 0.0005 if wmin > 0 else 0.001
    theHist = y.copy()
    if wmin > 0:
        theHist[0] = 0.0  # for the energy histograms

    bin_widths = np.diff(edges)
    s = np.sum(theHist * bin_widths)
    if s != 0:
        theHist = theHist / s

    am = np.argmax(theHist)
    MPV = 0.5 * (edges[am] + edges[am + 1])

    weights = theHist * bin_widths
    cumulative_weights_all = np.cumsum(weights)

    valid_centers = []
    valid_cumulative = []
    for i in range(len(edges) - 1):
        if edges[i + 1] < wmin or edges[i] > wmax:
            continue
        center = 0.5 * (edges[i + 1] + edges[i])
        valid_centers.append(center)
        valid_cumulative.append(cumulative_weights_all[i])

    def fallback_mean_stdev():
        centers_all = 0.5 * (edges[:-1] + edges[1:])
        mean = np.sum(centers_all * weights)
        std68 = np.sqrt(np.sum((centers_all - mean) ** 2 * weights))
        return std68, mean - std68, mean + std68, MPV

    if len(valid_centers) == 0:
        return fallback_mean_stdev()

    valid_centers = np.array(valid_centers)
    valid_cumulative = np.array(valid_cumulative)

    if NUMBA_AVAILABLE:
        width, low, high = _find_narrowest_interval(
            valid_centers, valid_cumulative, percentage, epsilon, wmin, wmax
        )
    else:
        width, low, high = _narrowest_interval_fallback(
            valid_centers, valid_cumulative, percentage, epsilon, wmin, wmax
        )

    if low == wmin and high == wmax:
        return fallback_mean_stdev()
    return 0.5 * (high - low), low, high, MPV


def sigma_rms(y, edges, wmin=0.7, wmax=1.2, **kwargs):
    yc = y.copy()
    MPV = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
    bin_widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    weights = yc * bin_widths
    norm = np.sum(weights)
    mean = np.sum(centers * weights) / norm
    RMS = np.sqrt(np.sum((centers - mean) ** 2 * weights) / norm)
    return RMS, mean - RMS, mean + RMS, MPV


def sigma_interquantile_range(y, edges, wmin=0.7, wmax=1.2, **kwargs):
    yc = y.copy()
    MPV = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
    s = np.sum(yc * np.diff(edges))
    if s != 0:
        yc = yc / s
    cumulative = np.cumsum(yc * np.diff(edges))
    low_idx = np.searchsorted(cumulative, 0.1585)
    high_idx = np.searchsorted(cumulative, 0.8415)
    low_idx = min(low_idx, len(edges) - 1)
    high_idx = min(high_idx, len(edges) - 1)
    low = edges[low_idx]
    high = edges[high_idx]
    return 0.5 * (high - low), low, high, MPV


def _double_crystal_ball(x, mu, sigma, alphaL, nL, alphaR, nR, norm):
    t = (x - mu) / sigma
    result = np.zeros_like(t)
    maskL = t < -alphaL
    result[maskL] = norm * (
        (nL / abs(alphaL)) ** nL
        * np.exp(-0.5 * alphaL**2)
        / (nL / abs(alphaL) - abs(alphaL) - t[maskL]) ** nL
    )
    maskC = (~maskL) & (t < alphaR)
    result[maskC] = norm * np.exp(-0.5 * t[maskC] ** 2)
    maskR = t >= alphaR
    result[maskR] = norm * (
        (nR / abs(alphaR)) ** nR
        * np.exp(-0.5 * alphaR**2)
        / (nR / abs(alphaR) - abs(alphaR) + t[maskR]) ** nR
    )
    return result


def sigma_dscb(y, edges, wmin=0.7, wmax=1.2, **kwargs):
    yc = y.copy()
    centers = 0.5 * (edges[1:] + edges[:-1])
    MPV_guess = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
    sigma_guess = np.std(np.repeat(centers, yc.astype(int))) if np.sum(yc) > 0 else 1.0
    p0 = [MPV_guess, sigma_guess, 1.5, 3.0, 1.5, 3.0, max(yc)]
    try:
        popt, _ = curve_fit(_double_crystal_ball, centers, yc, p0=p0, maxfev=10000)
    except RuntimeError:
        print("⚠️ DSCB fit failed.")
        return None
    mu, sigma, *_ = popt
    return sigma, mu - sigma, mu + sigma, mu


def sigma_gaussian_fit(y, edges, wmin=0.7, wmax=1.2, **kwargs):
    yc = y.copy()
    centers = 0.5 * (edges[1:] + edges[:-1])
    mask = (centers >= 0.75) & (centers <= 1.15)
    centers_fit = centers[mask]
    yc_fit = yc[mask]
    if len(centers_fit) < 3:
        print("⚠️ Not enough points to fit a Gaussian.")
        return None

    def gaussian(x, mu, sigma, norm):
        return norm * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    mean_guess = 0.5 * (centers_fit[np.argmax(yc_fit)] + centers_fit[np.argmax(yc_fit) + 1])
    sigma_guess = np.std(np.repeat(centers_fit, yc_fit.astype(int))) if np.sum(yc_fit) > 0 else 1.0
    p0 = [mean_guess, sigma_guess, max(yc_fit)]
    try:
        popt, _ = curve_fit(gaussian, centers_fit, yc_fit, p0=p0, maxfev=10000)
    except RuntimeError:
        print("⚠️ Gaussian fit failed.")
        return None
    mu, sigma, _ = popt
    return sigma, mu - sigma, mu + sigma, mu


SIGMA_METHODS = {
    "std68": sigma_std68,
    "RMS": sigma_rms,
    "interquantile_range": sigma_interquantile_range,
    "DSCB": sigma_dscb,
    "gaussian_fit": sigma_gaussian_fit,
}


def add_sigma_method(name, func):
    SIGMA_METHODS[name] = func


# --------------------------------------------------------------------------
# Energy/angle-dependence functional forms (fitted across energy/angle bins)
# --------------------------------------------------------------------------


def _model_three_param(E, a, b, c):
    return a / np.sqrt(E) + b + c / E


def _model_three_param_quadrature(E, a, b, c):
    return np.sqrt((a / np.sqrt(E)) ** 2 + b**2 + (c / E) ** 2)


def _model_two_param(E, a, b):
    return a / np.sqrt(E) + b


def _model_two_param_quadrature(E, a, b):
    return np.sqrt((a / np.sqrt(E)) ** 2 + b**2)


RESOLUTION_MODELS = {
    # sigma/E = A/sqrt(E) + B + C/E  (used for jet energy resolution)
    "three_param": dict(
        func=_model_three_param,
        p0=[0.5, 0.03, 0.1],
        bounds=([0.0, 0.0, 0.0], [np.inf, np.inf, np.inf]),
    ),
    "three_param_quadrature": dict(
        func=_model_three_param_quadrature,
        p0=[0.5, 0.03, 0.1],
        bounds=([0.0, 0.0, 0.0], [np.inf, np.inf, np.inf]),
    ),
    # sigma = A/sqrt(E) + B  (used for angular resolution)
    "two_param": dict(
        func=_model_two_param,
        p0=[0.5, 0.03],
        bounds=([0.0, 0.0], [np.inf, np.inf]),
    ),
    "two_param_quadrature": dict(
        func=_model_two_param_quadrature,
        p0=[0.5, 0.03],
        bounds=([0.0, 0.0], [np.inf, np.inf]),
    ),
}


def add_resolution_model(name, func, p0, bounds):
    RESOLUTION_MODELS[name] = dict(func=func, p0=p0, bounds=bounds)


def fit_resolution_model(mid_points, values, model="three_param", min_E=0.0,
                          bounds_override=None, n_curve_points=100):
    """Fit a registered energy/angle-dependence model.

    Returns (xs, ys, popt, pcov) where xs/ys trace the fitted curve across
    the range of mid_points, and popt/pcov are the fit parameters/covariance
    (popt sign-corrected via abs(), matching the original behaviour).
    """
    spec = RESOLUTION_MODELS[model]
    mid_points = np.asarray(mid_points, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = mid_points >= min_E
    mid_points = mid_points[mask]
    values = values[mask]
    bounds = bounds_override if bounds_override is not None else spec["bounds"]
    popt, pcov = curve_fit(
        spec["func"], mid_points, values, p0=spec["p0"], maxfev=10000, bounds=bounds
    )
    xs = np.linspace(mid_points.min(), mid_points.max(), n_curve_points)
    ys = spec["func"](xs, *popt)
    return xs, ys, np.abs(popt), pcov


def downsample_for_dashboard(y, edges, max_points=150):
    """Rebin (y, edges) for storage in the (small) dashboard data file.

    Purely for visual display in the dashboard - not used for any
    statistical computation, which always happens on the full-resolution
    histograms. Generic over any binned histogram, so it's shared between
    the resolution and Higgs-mass dashboard pickles.
    """
    y = np.asarray(y)
    edges = np.asarray(edges)
    n = len(y)
    if n <= max_points:
        return y, edges
    factor = int(np.ceil(n / max_points))
    n_new = n // factor
    y_ds = np.array([y[i * factor : (i + 1) * factor].sum() for i in range(n_new)])
    edges_ds = np.array([edges[i * factor] for i in range(n_new)] + [edges[n_new * factor]])
    return y_ds, edges_ds


# --------------------------------------------------------------------------
# Peak-shape fit models (used for the Higgs mass peak; pluggable so more
# fitting methods can be added later without touching mass_plots.py).
# --------------------------------------------------------------------------


def _gaussian_peak(x, mu, sigma, norm):
    return norm * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


PEAK_FIT_MODELS = {
    "gaussian": dict(func=_gaussian_peak, n_params=3),
    "dscb": dict(func=_double_crystal_ball, n_params=7),
}


def add_peak_fit_model(name, func, n_params):
    PEAK_FIT_MODELS[name] = dict(func=func, n_params=n_params)


def fit_peak(x_vals, y_vals, model="gaussian", window=None, n_curve_points=200):
    """Fit a peak-shape model to a 1D histogram (e.g. the Higgs mass peak).

    Restricts the fit to `window` (mu0 +/- 30, by default, where mu0 is the
    argmax of y_vals) so a long tail doesn't bias the fit. Returns
    (xs, ys, popt) tracing the fitted curve, or None on failure/too few points.
    """
    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    if len(x) == 0 or np.sum(y) <= 0:
        return None

    mu0 = x[np.argmax(y)]
    if window is None:
        window = (mu0 - 30, mu0 + 30)
    mask = (x >= window[0]) & (x <= window[1])
    x_fit, y_fit = x[mask], y[mask]
    if len(x_fit) < 5:
        return None

    spec = PEAK_FIT_MODELS[model]
    sigma0 = max(float(np.sqrt(np.average((x_fit - mu0) ** 2, weights=np.clip(y_fit, 0, None)))), 1e-3)
    norm0 = float(np.max(y_fit))

    if model == "gaussian":
        p0 = [mu0, sigma0, norm0]
        bounds = ([x_fit.min(), 1e-3, 0.0], [x_fit.max(), (x_fit.max() - x_fit.min()), np.inf])
    elif model == "dscb":
        p0 = [mu0, sigma0, 1.5, 3.0, 1.5, 3.0, norm0]
        bounds = (-np.inf, np.inf)
    else:
        raise ValueError(f"Unknown peak fit model: {model}")

    try:
        popt, _ = curve_fit(spec["func"], x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
    except RuntimeError:
        return None

    xs = np.linspace(x_fit.min(), x_fit.max(), n_curve_points)
    ys = spec["func"](xs, *popt)
    return xs, ys, popt
