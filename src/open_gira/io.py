"""
Common functionality for reading and writing to disk.
"""

import functools
from glob import glob
import logging
import json
from os.path import splitext, basename, join
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from tqdm import tqdm
import xarray as xr

from open_gira.utils import natural_sort


WGS84_EPSG = 4326

# lines beginning with this character will be ignored by pandas
COMMENT_PREFIX: str = "#"


def bit_pack_dataset_encoding(ds: xr.Dataset, n_bits: int = 16) -> dict:
    """
    Given a Dataset, return the dict to pack its variables into an integer type
    for serialisation as netCDF. Should provide space savings over using a
    single or double floating point type.

    Args:
        ds: Dataset to determine encoding parameters for.
        n_bits: Size of integer type to use as storage on disk.

    Returns:
        Encoding dictionary for use with xr.Dataset.to_netcdf.
    """

    if not isinstance(ds, xr.Dataset):
        raise ValueError(f"ds must be an xarray Dataset, but is instead {type(ds)}")

    per_variable: list[dict] = [
        bit_pack_dataarray_encoding(ds[var_name], n_bits) for var_name in ds.data_vars
    ]

    # merge the dictionaries from the DataArrays
    return functools.reduce(lambda a, b: {**a, **b}, per_variable)


def bit_pack_dataarray_encoding(da: xr.DataArray, n_bits: int = 16) -> dict:
    """
    Given a DataArray, return the dict to pack it into an integer type for
    serialisation as netCDF. Should provide space savings over using a single or
    double floating point type.

    Args:
        da: DataArray to determine encoding parameters for.
        n_bits: Size of integer type to use as storage on disk.

    Returns:
        Encoding dictionary for use with xr.DataArray.to_netcdf.
    """

    if not isinstance(da, xr.DataArray):
        raise ValueError(f"da must be an xarray DataArray, but is instead {type(da)}")

    def defaults() -> Tuple[int, int, int]:
        return 1, 0, -1

    if da.size == 0:
        # no data
        scale_factor, add_offset, fill_value = defaults()
    elif np.count_nonzero(np.isnan(da.values)) == da.size:
        # everything is NaN!
        scale_factor, add_offset, fill_value = defaults()
    else:
        scale_factor, add_offset, fill_value = netcdf_packing_parameters(
            da.min().item(),
            da.max().item(),
            n_bits
        )

    return {
        da.name: {
            'dtype': f'int{n_bits:d}',
            'scale_factor': scale_factor,
            'add_offset': add_offset,
            '_FillValue': fill_value
        }
    }


def netcdf_packing_parameters(minimum: float, maximum: float, n_bits: int) -> Tuple[float, float]:
    """
    Given (floating point) data within a certain range, find the best scale
    factor and offset to use to pack as signed integer values, using most of
    the available bits:

        deserialised = scale_factor * serialised + add_offset
        serialised = (deserialised - add_offset) / scale_factor

    The fill_value is a sentinel value required to successfully round-trip NaN
    values from float (in memory) to int (on disk) to float (in memory).

    The precision of the round-tripped data is the scale_factor, i.e. for a
    scale_factor of 0.01, deserialised figures will differ by increments of
    0.01. If this is approximately equal to (or larger than) features of
    interest in your data, you will lose those features.

    N.B. Whatever reads this data from disk must read and employ the
    scale_factor and add_offset metadata! netCDF4-python and xarray do, ncdump
    doesn't.

    See here for more information:
    https://docs.xarray.dev/en/stable/user-guide/io.html#scaling-and-type-conversions

    Args:
        maximum: Maximum value in data to serialise
        minimum: Minimum value in data to serialise
        n_bits: Number of available bits

    Returns:
        scale_factor: Used to (de)serialise: decoded = scale_factor * encoded + add_offset
        add_offset: Used to (de)serialise: decoded = scale_factor * encoded + add_offset
        fill_value: Sentinel value in integer space used for storing NaN
    """

    # need a finite range
    assert np.isfinite(minimum)
    assert np.isfinite(maximum)
    assert maximum - minimum >= 0

    # _FillValue used to representing NaN as serialised integer
    # we have kept room at the ends of the integer bit space to avoid a collision
    fill_value = -2 ** (n_bits - 1)

    # if there is no variance in the data, return unscaled
    if minimum == maximum:
        # use a fill value that doesn't collide with the data
        if minimum == maximum == fill_value:
            return 1, 0, -1 * fill_value
        else:
            return 1, 0, fill_value

    # use this fraction of the available `n_bits` space
    # this leaves room at either end (-2**n_bits or 2**n_bits) for _FillValue
    occupancy_fraction = 0.99

    # factor to rescale data to (almost) the available 2 ** n_bits serialised range
    scale_factor = (maximum - minimum) / (2 ** (occupancy_fraction * n_bits))

    # offset to center serialised data about 0
    add_offset = minimum + 2 ** (occupancy_fraction * n_bits - 1) * scale_factor

    return scale_factor, add_offset, fill_value


@functools.lru_cache(maxsize=128)
def cached_json_file_read(file_path: str):
    """
    Read a JSON file and return its parsed contents. Cached on function argument.

    Args:
        file_path: Path to JSON file to read.

    Returns:
        Parsed contents of `file_path`.
    """

    with open(file_path, "r") as fp:
        return json.load(fp)


def concat_geoparquet(paths: list[str]) -> gpd.GeoDataFrame:
    """
    Sort list of paths, create GeoDataFrames, concatenate into a single
    GeoDataFrame, reset the index and return.
    """

    dataframes: list[gpd.GeoDataFrame] = []

    # sort paths for reproducibility
    for i, path in tqdm(enumerate(natural_sort(paths))):

        try:
            gdf = gpd.read_parquet(path)

        except ValueError as error:
            if NO_GEOM_ERROR_MSG in str(error):
                # if the input parquet file does not contain a geometry column,
                # geopandas will raise a ValueError rather than try to procede. we
                # catch that here, but check the error message - to be more
                # specific than catching and suppressing any ValueError

                # use an empty geodataframe to append instead
                gdf = gpd.GeoDataFrame([])

        dataframes.append(gdf)

    logging.info("Joining GeoDataFrames")

    # pandas concat of iterable containing GeoDataFrames will return a GeoDataFrame
    concatenated = pd.concat(dataframes)

    return concatenated.reset_index(drop=True)


def write_empty_frames(edges_path: str, nodes_path: Optional[str] = None) -> None:
    """
    If we don't have sufficient / good enough input data, write out empty output.

    N.B. Output files must exist for snakemake's sake.
    """

    # write with a CRS, makes it easier to concatenate dataframes later
    empty_gdf = gpd.GeoDataFrame({"geometry": []}, crs=pyproj.CRS.from_user_input(WGS84_EPSG))
    empty_gdf.to_parquet(edges_path)

    # some parts of the workflow only consider edges, not nodes
    # when not passed a nodes_path, do not attempt to write
    if nodes_path:
        empty_gdf.to_parquet(nodes_path)

    return


def read_damage_curves(damage_curves_dir: str, hazard_type: str, asset_types: set) -> dict[str, pd.DataFrame]:
    """
    Load damage curves from CSVs on disk

    Expected to reside in following structure:
    <damage_curves_dir>/<hazard_type>/<asset_type>.csv

    Damage curve files may have comments, these are lines starting with COMMENT_PREFIX

    Args:
        damage_curves_dir (str): Path to folder containing hazards
        hazard_type (str): Name of hazard folder containing asset specific curves
        asset_types (set): Asset types we require damage curves for

    Returns (dict[str, pd.DataFrame):
        Mapping from asset_type to respective damage curve
    """

    # fetch damage curves for relevant hazard type
    damage_curve_paths = glob(join(damage_curves_dir, hazard_type, "*.csv"))

    damage_curves: dict[str, pd.DataFrame] = {
        # curves expected to be named as a value of Asset class, e.g. RoadAssets.BRIDGE -> road_bridge.csv
        # dict is asset_type: dataframe with hazard intensity [0, inf] and damage fraction [0, 1]
        splitext(basename(path))[0]: pd.read_csv(path, comment=COMMENT_PREFIX) for path in damage_curve_paths
    }

    for asset_type, damage_curve in damage_curves.items():
        assert len(damage_curve.columns) == 2
        # check hazard intensity and damage fraction values are 0 or positive real
        assert ((damage_curve >= 0).all()).all()
        # check damage fraction is less than or equal to 1
        assert (damage_curve.iloc[:, 1] <= 1).all()

    if not set(damage_curves.keys()).issuperset(asset_types):
        raise RuntimeError(f"requested {asset_types=} not all found: {damage_curves.keys()=}")

    return damage_curves


def read_rehab_costs(path: str) -> pd.DataFrame:
    """
    Read a rehabilitation costs CSV table and check its contents. Expects two columns:
        1) `asset_type`, column containing categorical strings, as enumerated in open_gira.assets.
        2) Some sort of monetary cost in US Dollars. What this cost refers to exactly may depend on the context.

    Arguments:
        path: Path to rehabilitation cost CSV file.

    Returns:
        Rehabilitation costs
    """
    costs = pd.read_csv(path, comment=COMMENT_PREFIX)

    assert len(costs) > 0
    assert len(costs.columns) == 2

    # check asset_type
    assert 'asset_type' == costs.columns[0]
    assert costs.asset_type.dtype == object

    # check costs
    assert costs.iloc[:, 1].dtype == float
    assert (costs.iloc[:, 1] >= 0).all()

    return costs