"""
Common functionality for reading and writing to disk.
"""

from glob import glob
from os.path import splitext, basename, join
from typing import Optional

import geopandas as gpd
import pandas as pd
import pyproj


WGS84_EPSG = 4326


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

    # lines beginning with this character will be ignored by pandas
    COMMENT_PREFIX: str = "#"

    # fetch damage curves for relevant hazard type
    damage_curve_paths = glob(join(damage_curves_dir, hazard_type, "*.csv"))

    damage_curves: dict[str, pd.DataFrame] = {
        # curves expected to be named as a value of Asset class, e.g. RoadAssets.BRIDGE -> road_bridge.csv
        # dict is asset_type: dataframe with hazard intensity [0, inf] and damage fraction [0, 1]
        splitext(basename(path))[0]: pd.read_csv(path, comment=COMMENT_PREFIX) for path in damage_curve_paths
    }

    for asset_type, damage_curve in damage_curves.items():
        # check hazard intensity and damage fraction values are 0 or positive real
        assert ((damage_curve >= 0).all()).all()
        # check damage fraction is less than or equal to 1
        assert (damage_curve.iloc[:, 1] <= 1).all()

    if not set(damage_curves.keys()).issuperset(asset_types):
        raise RuntimeError(f"requested {asset_types=} not all found: {damage_curves.keys()=}")

    return damage_curves
