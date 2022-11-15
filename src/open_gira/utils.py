"""
Reusable functionality that doesn't clearly fit anywhere else.
"""


import re
from typing import Iterable

import pandas as pd


def natural_sort(to_sort: Iterable) -> list:
    """
    Natural sort iterables of strings by splitting on digits and sorting by the
    resulting list.
    """

    return sorted(
        to_sort,
        key=lambda entry: [
            int(fragment) if fragment.isdigit() else fragment.lower()
            for fragment in re.split(r"(\d+)", entry)
        ]
    )


def str_to_bool(series: pd.Series) -> pd.Series:
    """
    Make a stab at turning a series of strings into a boolean series.

    Args:
        series (pd.Series): Input series

    Returns:
        pd.Series: Boolean series of whether or not strings are truthy (by our logic).
    """

    FALSE_VALUES = {'n', 'no', 'false', 'f', ''}

    def str_parser(s: str) -> bool:
        """
        If a string is in the set below, return False, otherwise, return True.
        """

        if not isinstance(s, str):
            raise ValueError(f"{s=} has {type(s)=}, but should be str")

        return s.lower() not in FALSE_VALUES

    # set our null values to the empty string
    new_series = series.copy(deep=True)
    new_series.loc[series.isnull()] = ''

    return new_series.apply(str_parser)
