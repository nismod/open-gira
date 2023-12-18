"""
Tools for manipulating administrative region data.
"""

import logging

import pandas as pd


def merge_gadm_admin_levels(preference: pd.DataFrame, alternative: pd.DataFrame) -> pd.DataFrame:
    """
    Geospatial data is often aggregated at admin 'levels', as exemplified by the
    GADM project. These integer levels are 0 for nation states, 1 for the
    regions of the state, 2 for sub-regions and so on. This function is designed
    to take two different levels, e.g. 2 and 1 as `preference` and
    `alternative`, respectively and combine them. Where a country has
    representation in `preference`, these rows will be returned by default. If
    there is a country in `alternative` that is not in `preference`, this
    country's rows from `alternative` will be merged into the returned result.

    Args:
        preference: Table with rows to return by default. Must contain an ISO_A3
            column.
        alternative: Table with rows to choose from if `preference` is missing a
            country. Must contain an ISO_A3 column.

    Returns:
        Table of rows representing the union of countries in the passed tables,
        with data from `preference` by default, but `alternative` if `preference`
        data is not available.
    """

    substitute_countries = set(alternative.ISO_A3) - set(preference.ISO_A3)
    logging.info(f"Gap filling with: {substitute_countries}")

    substitute_regions: pd.DataFrame = alternative[alternative["ISO_A3"].isin(substitute_countries)]

    merged = pd.concat([preference, substitute_regions])

    return merged.sort_values("ISO_A3")