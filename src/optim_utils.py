# -*- coding: utf-8 -*-
"""
Shared scalar root-finding helper used by the pipe and heat-exchanger models.

Several components solve a 1D energy-balance equation every timestep by
minimizing a squared residual with `scipy.optimize.minimize_scalar` (bounded
golden-section search, typically 30-100 function evaluations). Since the
residual is a true root-finding problem (not a genuine minimization), solving
it with `scipy.optimize.brentq` on the raw (non-squared) residual converges
in far fewer evaluations, and warm-starting the bracket around the previous
timestep's solution reduces that further still.
"""
import scipy.optimize as opt


def robust_root_scalar(residual_fn, lo, hi, prev_guess=None, narrow=0.5, xtol=1e-6):
    '''
    Solve residual_fn(x) == 0 within [lo, hi] using brentq.

    Tries a narrow bracket around prev_guess first (warm start from the
    previous timestep's solution), falling back to the full [lo, hi]
    bracket, and finally to a bounded minimize_scalar of the squared
    residual if no sign change is found anywhere (safety net that preserves
    the old, more permissive behaviour for edge cases where the residual
    does not bracket a root, e.g. right at the operating limits).

    Evaluating residual_fn exactly at lo/hi can occasionally hit a
    thermodynamic singularity (e.g. CoolProp raising right at a saturation
    boundary) that the old bounded minimize_scalar never actually probed
    (it only evaluates interior points) - any such exception is treated as
    "no usable bracket here" and we fall through to the next strategy.
    '''
    if lo > hi:
        lo, hi = hi, lo

    def safe_eval(x):
        try:
            return residual_fn(x)
        except Exception:
            return None

    if prev_guess is not None and lo <= prev_guess <= hi:
        lo_n, hi_n = max(lo, prev_guess - narrow), min(hi, prev_guess + narrow)
        if lo_n < hi_n:
            f_lo, f_hi = safe_eval(lo_n), safe_eval(hi_n)
            if f_lo is not None and f_hi is not None:
                if f_lo == 0.0:
                    return lo_n
                if f_hi == 0.0:
                    return hi_n
                if f_lo * f_hi < 0.0:
                    try:
                        return opt.brentq(residual_fn, lo_n, hi_n, xtol=xtol)
                    except Exception:
                        pass

    f_lo, f_hi = safe_eval(lo), safe_eval(hi)
    if f_lo is not None and f_hi is not None:
        if f_lo == 0.0:
            return lo
        if f_hi == 0.0:
            return hi
        if f_lo * f_hi < 0.0:
            try:
                return opt.brentq(residual_fn, lo, hi, xtol=xtol)
            except Exception:
                pass

    def squared(x):
        return residual_fn(x) ** 2
    rez = opt.minimize_scalar(squared, bounds=(lo, hi), method='bounded')
    return rez.x
