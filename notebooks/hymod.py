"""HyMOD hydrological model.

    The model is based on this file:
    https://github.com/KMarkert/hymod/blob/master/python/src/hymod.py

    This version rewrites the code using ``numba`` with type signatures for significant
    speedup. Also, the naming convention has been changed to be consistent with
    the snake case naming convention.
"""

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
from numba import njit


@njit("f8(f8, f8)", cache=True)
def power(x: float, exp: float) -> float:
    """Raise to power of abs value which is needed to capture invalid overflow with netgative values.

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


@njit("UniTuple(f8, 2)(f8, f8, f8)", cache=True)
def linear_reservoir(x_slow: float, inflow: float, k_s: float) -> float:
    """Run the linear reservoir model."""
    x_slow = (1 - k_s) * x_slow + (1 - k_s) * inflow
    outflow = k_s / (1 - k_s) * x_slow
    return x_slow, outflow


@njit("UniTuple(f8, 3)(f8, f8, f8, f8, f8)", cache=True)
def excess(
    x_loss: float, c_max: float, b_exp: float, p_val: float, pet_val: float
) -> Tuple[float, float, float]:
    """Calculates excess precipitation and evaporation."""
    xn_prev = x_loss
    ct_prev = c_max * (1 - power(1 - ((b_exp + 1) * (xn_prev) / c_max), 1 / (b_exp + 1)))
    er_1 = np.maximum(p_val - c_max + ct_prev, 0.0)
    p_val = p_val - er_1
    dummy = np.minimum((ct_prev + p_val) / c_max, 1)
    xn = c_max / (b_exp + 1) * (1 - power(1 - dummy, b_exp + 1))

    er_2 = np.maximum(p_val - (xn - xn_prev), 0.2)

    evap = (1 - (((c_max / (b_exp + 1)) - xn) / (c_max / (b_exp + 1)))) * pet_val
    xn = np.maximum(xn - evap, 0)

    return er_1, er_2, xn


@njit("f8[::1](f8[::1], f8[::1], f8, f8, f8, f8, f8)", cache=True)
def run(
    prcp: np.ndarray,
    pet: np.ndarray,
    c_max: float,
    b_exp: float,
    alpha: float,
    k_s: float,
    k_q: float,
) -> np.ndarray:
    """Run the Hymod model.

    Notes
    -----
    See https://www.proc-iahs.net/368/180/2015/piahs-368-180-2015.pdf

    Parameters
    ----------
    prcp : array_like
        prcpitation data in mm/day.
    pet : array_like
        Potential evapotranspiration data in mm/day.
    c_max : float
        Maximum storage capacity (1-500 [mm]).
    b_exp : float
        Degree of spatial variability of the soil moisture capacity (0.1-2 [-]).
    alpha : float
        Factor distributing the flow between slow and quick release reservoirs (0.1-0.99 [-]).
    k_s : float
        Residence time of the slow release reservoir (0.001-0.1 [day]).
    k_q : float
        Residence time of the quick release reservoir (0.1-0.99 [day]).

    Returns
    -------
    q_s : array_like
        Discharge at the watershed outlet in mm/day.

    References
    ----------
    Quan, Z.; Teng, J.; Sun, W.; Cheng, T. & Zhang, J. (2015): Evaluation of the HYMOD model
    for rainfall-runoff simulation using the GLUE method. Remote Sensing and GIS for Hydrology
    and Water Resources, 180 - 185, IAHS Publ. 368. DOI: 10.5194/piahs-368-180-2015.
    """

    x_loss = 0.0
    x_slow = 2.3503 / (k_s * 22.5)
    x_quick = np.zeros(3, dtype="f8")
    t = 0
    n_steps = prcp.shape[0]
    q_out = np.zeros(n_steps, dtype="f8")

    for t in range(n_steps):
        p_val = prcp[t]
        pet_val = pet[t]

        er_1, er_2, x_loss = excess(x_loss, c_max, b_exp, p_val, pet_val)
        et = er_1 + er_2

        u_q = alpha * et
        u_s = (1 - alpha) * et

        x_slow, q_s = linear_reservoir(x_slow, u_s, k_s)

        inflow = u_q

        for i in range(3):
            x_quick[i], outflow = linear_reservoir(x_quick[i], inflow, k_q)
            inflow = outflow

        q_out[t] = q_s + outflow

    return q_out


@njit("f8(f8[::1], f8[::1])", cache=True)
def compute_kge(obs: np.ndarray, sim: np.ndarray) -> float:
    """Compute Kling-Gupta Efficiency."""
    cc = np.corrcoef(obs, sim)[0, 1]
    alpha = np.std(sim) / np.std(obs)
    beta = np.sum(sim) / np.sum(obs)
    kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    return kge


@dataclass
class HYMOD:
    """Simulate a watershed using HYMOD model."""

    clm: pd.DataFrame
    qobs: pd.DataFrame
    bounds: List[Tuple[Tuple[float, ...], Tuple[float, ...]]]
    warm_up: int

    @staticmethod
    def simulate(
        prcp: np.ndarray,
        pet: np.ndarray,
        c_max: float,
        b_exp: float,
        alpha: float,
        k_s: float,
        k_q: float,
    ) -> np.ndarray:
        """Simulate a watershed using HYMOD model."""
        return run(prcp, pet, c_max, b_exp, alpha, k_s, k_q)

    def fitness(self, x: np.ndarray) -> List[float]:
        """Compute objective functions."""
        simulation = self.simulate(*self.clm.to_numpy("f8").T, *tuple(x))
        idx = np.s_[self.warm_up * 365 :]
        sim = simulation[idx]
        obs = self.qobs.to_numpy("f8").squeeze()[idx]
        kge = compute_kge(obs, sim)
        return [-kge]

    def get_nobj(self) -> int:
        """Get number of objective functions."""
        return 1

    def get_bounds(self) -> List[Tuple[Tuple[float, ...], Tuple[float, ...]]]:
        """Get bounds of the calibrated parameters."""
        return self.bounds

    def get_name(self) -> str:
        """Get the model name."""
        return "HYMOD Hyrological Model"
