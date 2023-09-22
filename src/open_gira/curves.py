"""
Functions describing mathematical curves.
"""


import numpy as np


def logistic_min(x: float | np.ndarray, L: float, m: float, k: float, x_0: float) -> float | np.ndarray:
    """
    Logistic function with a minimum value, m.

    Args:
        x: Input values
        L: Maximum output value
        m: Minimum output value
        k: Steepness parameter
        x_0: Location of sigmoid centre in x

    Returns:
        Output values
    """

    return m + (L - m) / (1 + np.exp(-k * (x - x_0)))
