"""Analysis functions for common FLASH problems."""

import numpy as np
import scipy.optimize


def fit_power_law(t, r, **kwargs):
    """Given a shock with radius r over time t, returns the parameters of best fit
    to the following power-law:

    .. math:: R = R_0 t^\\alpha

    Parameters
    ----------
    t : iterable
        Time values.
    r : iterable
        Corresponding radius values.
    kwargs : optional
        Other parameters to pass down to scipy.optimize.curve_fit()

    Returns
    -------
    r_o : float
        Radius scaling factor.
    alpha : float
        The power-law index.

    """
    power_law = lambda s, r0, a: r0 * (s**a)
    popt, pcov = scipy.optimize.curve_fit(power_law, t, r, p0=(r[0], 0.5), **kwargs)
    return popt
    

def shock_detect(r, v, threshold=1e-6, min_threshold=1e-36):
    """Given radius (r) and data (v) arrays, finds the outer most shock.
    Currently this assumes that the radius is montonically increasing.

    Parameters
    ----------
    r : ndarray
        Positional data for values.
    v : ndarray
        Data values which determine the presnece of a shock.  This is 
        typically done with density or pressure fields.
    threshold : float, optional
        The value above which gradients must exceed to be considered a
        shock.  Required for ignoring low-level noise.
    min_threshold : float, optional
        If a shock is not found at a given threshold level, the threshold
        value is reduced by an order of magnitude until a shock is found.
        This continues until a minimum threshold is reached.

    Returns
    -------
    shock_r : r.dtype
        Shock position.
    shock_v : v.dtype
        Shock peak value.
    shock_i : int
        Shock index in r and v arrays

    """
    dvdr = np.gradient(v) / np.gradient(r)
    abs_v = np.abs(v[:-1])
    sign_dvdr = np.sign(dvdr)
    mask = np.array([])
    while 0 == mask.sum() and min_threshold <= threshold:
        mask = (threshold < np.abs(v[:-1]))
        mask = mask & (sign_dvdr[1:] != sign_dvdr[:-1])
        threshold *= 0.1
    if threshold < min_threshold:
        raise ValueError("minimum threshold reached, shock not found.")
    shock_i = np.where(mask)[0][-1]
    shock_r = r[shock_i]
    shock_v = v[shock_i]
    return shock_r, shock_v, shock_i
