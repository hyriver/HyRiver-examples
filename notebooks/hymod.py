"""HyMOD hydrological model.

The model is based on this file:
https://github.com/KMarkert/hymod/blob/master/python/src/hymod.py

This version rewrites the code using ``numba`` with type signatures for significant
speed up. Also, the naming convention has been changed to be consistent with
the snake case naming convention.
"""
import functools
import warnings
from typing import Any, Dict, List, Tuple

import numpy as np
import numpy.typing as npt
import pandas as pd

try:
    from numba import njit

    ngjit = functools.partial(njit, cache=True, nogil=True)
except ImportError:
    warnings.warn("Numba not installed. Using slow pure python version.", UserWarning)

    def ngjit(_):
        def decorator_njit(func):
            @functools.wraps(func)
            def wrapper_decorator(*args, **kwargs):
                return func(*args, **kwargs)

            return wrapper_decorator

        return decorator_njit


ZERO = np.float64(0)


@ngjit("f8(f8, f8)")
def power(x: np.float64, exp: np.float64) -> np.float64:
    """Raise to power of abs value needed to capture invalid overflow with netgative values.

    Parameters
    ----------
    x : float
        Input array.
    exp : float
        Exponent.

    Returns
    -------
    result : array_like
        Result of the power operation.
    """
    return np.power(np.abs(x), exp)


@ngjit("UniTuple(f8, 2)(f8, f8, f8)")
def linear_reservoir(
    x_slow: np.float64, inflow: np.float64, k_s: np.float64
) -> Tuple[np.float64, np.float64]:
    """Run the linear reservoir model."""
    xn_slow = (1 - k_s) * x_slow + (1 - k_s) * inflow
    outflow = k_s / (1 - k_s) * xn_slow
    return xn_slow, outflow


@ngjit("UniTuple(f8, 2)(f8, f8, f8, f8, f8)")
def excess(
    prcp_t: np.float64,
    pet_t: np.float64,
    x_loss: np.float64,
    c_max: np.float64,
    b_exp: np.float64,
) -> Tuple[np.float64, np.float64]:
    """Calculates excess precipitation and evaporation."""
    ct_prev = c_max * (1 - power(1 - ((b_exp + 1) * (x_loss) / c_max), 1 / (b_exp + 1)))
    er_1 = np.maximum(prcp_t - c_max + ct_prev, ZERO)
    s_1 = prcp_t - er_1
    dummy = np.minimum((ct_prev + s_1) / c_max, 1)
    xn = c_max / (b_exp + 1) * (1 - power(1 - dummy, b_exp + 1))

    er_2 = np.maximum(s_1 - (xn - x_loss), 0.2)

    evap = (1 - (((c_max / (b_exp + 1)) - xn) / (c_max / (b_exp + 1)))) * pet_t
    xn = np.maximum(xn - evap, 0)

    return er_1 + er_2, xn


@ngjit("f8[::1](f8[::1], f8[::1], f8, f8, f8, f8, f8)")
def run(
    prcp: npt.NDArray[np.float64],
    pet: npt.NDArray[np.float64],
    c_max: np.float64,
    b_exp: np.float64,
    alpha: np.float64,
    k_s: np.float64,
    k_q: np.float64,
) -> npt.NDArray[np.float64]:
    """Run the Hymod [1] model.

    Parameters
    ----------
    prcp : array_like
        prcpitation data in mm/day.
    pet : array_like
        Potential evapotranspiration data in mm/day.
    c_max : float
        Maximum storage capacity (1-1500 [mm]).
    b_exp : float
        Degree of spatial variability of the soil moisture capacity (0-2 [-]).
    alpha : float
        Factor distributing the flow between slow and quick release reservoirs (0.2-0.99 [-]).
    k_s : float
        Residence time of the slow release reservoir (0.01-0.5 [day]).
    k_q : float
        Residence time of the quick release reservoir (0.5-1.2 [day]).

    Returns
    -------
    array_like
        Discharge at the watershed outlet in mm/day.

    References
    ----------
    [1] Quan, Z.; Teng, J.; Sun, W.; Cheng, T. & Zhang, J. (2015): Evaluation of the HYMOD model
    for rainfall-runoff simulation using the GLUE method. Remote Sensing and GIS for Hydrology
    and Water Resources, 180 - 185, IAHS Publ. 368. DOI: 10.5194/piahs-368-180-2015.
    """
    x_loss = ZERO
    x_slow = ZERO
    x_quick = np.zeros(3, dtype="f8")
    n_steps = prcp.shape[0]
    q_out = np.zeros(n_steps, dtype="f8")

    for t in range(n_steps):
        et, x_loss = excess(prcp[t], pet[t], x_loss, c_max, b_exp)

        u_q = alpha * et
        u_s = (1 - alpha) * et

        x_slow, q_s = linear_reservoir(x_slow, u_s, k_s)

        inflow = u_q

        for i in range(3):
            x_quick[i], outflow = linear_reservoir(x_quick[i], inflow, k_q)
            inflow = outflow

        q_out[t] = q_s + outflow

    return q_out


@ngjit("f8(f8[::1], f8[::1])")
def compute_kge(sim: npt.NDArray[np.float64], obs: npt.NDArray[np.float64]) -> np.float64:
    """Compute Kling-Gupta Efficiency."""
    cc = np.corrcoef(obs, sim)[0, 1]
    alpha = np.std(sim) / np.std(obs)
    beta = np.sum(sim) / np.sum(obs)
    kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    return kge


class HYMOD:
    """Simulate a watershed using HYMOD model."""

    def __init__(self, clm: pd.DataFrame, qobs: pd.Series, cmax: float, warm_up: int) -> None:
        """Initialize the model.

        Parameters
        ----------
        clm : pandas.DataFrame
            Climate data with two columns: ``prcp`` and ``pet``.
        qobs : pandas.Series
            Streamflow observations.
        warm_up : int
            Number of warm-up years.
        """
        self.prcp = clm.prcp.to_numpy("f8")
        self.pet = clm.pet.to_numpy("f8")
        self.qobs = qobs.to_numpy("f8")
        self.cmax = np.float64(cmax)
        self.cal_idx = np.s_[warm_up * 365 :]

    def simulate(self, x: npt.NDArray[np.float64]) -> np.float64:
        """Compute objective functions."""
        b_exp, alpha, k_s, k_q = x
        return run(self.prcp, self.pet, self.cmax, b_exp, alpha, k_s, k_q)

    def fit(self, x: npt.NDArray[np.float64]) -> np.float64:
        """Compute objective functions."""
        b_exp, alpha, k_s, k_q = x
        qsim = run(self.prcp, self.pet, self.cmax, b_exp, alpha, k_s, k_q)
        kge = compute_kge(qsim[self.cal_idx], self.qobs[self.cal_idx])
        return -kge
