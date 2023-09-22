"""
For a given country and storm set, combine all the storm disruption estimates and
aggregate them to some administrative level.

Write out expected annual disruption (population of region affected by loss of
service averaged over all events).
"""


import logging
import sys

import geopandas as gpd
import pandas as pd
import xarray as xr


def disruption_by_target(disruption_path: str) -> tuple[str, pd.DataFrame]:
    """
    Given a disruption file, read netCDF and transform to pandas dataframe.

    Args:
        disruption_path: Path to a disruption netCDF containing `customers_affected`
            variable on `event_id` (singleton), `threshold` and `target`
            coordinates.

    Returns:
        Event name
        Disruption dataframe with target rows and threshold columns for given
            storm. Values are numbers of people impacted by power losses.
    """

    disruption = xr.open_dataset(disruption_path)

    event_id: str = disruption.event_id.item()
    logging.info(event_id)

    df: pd.DataFrame = disruption.customers_affected.to_pandas().T

    # N.B. arrow specification requires string column names (but threshold column names are float)
    df.columns = df.columns.astype(str)

    return event_id, df


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # load aggregation regions for level and country in question
    logging.info("Loading regions")
    admin: gpd.GeoDataFrame = gpd.read_parquet(snakemake.input.admin_areas)
    admin_level = int(snakemake.wildcards.ADMIN_SLUG.split("-")[-1])
    regions: gpd.GeoDataFrame = \
        admin[admin.GID_0 == snakemake.wildcards.COUNTRY_ISO_A3][[f"NAME_{admin_level}", f"GID_{admin_level}", "geometry"]]

    # load tracks (we will lookup storm dates from here)
    logging.info("Loading tracks")
    tracks: pd.DataFrame = pd.read_parquet(snakemake.input.tracks, columns=["track_id", "year"])
    track_year: pd.DataFrame = tracks.drop_duplicates("track_id").set_index("track_id")

    # load country targets file
    logging.info("Loading targets")
    targets: gpd.GeoDataFrame = gpd.read_parquet(snakemake.input.targets).set_index("id", drop=True)

    n_proc: int = snakemake.threads
    logging.info(f"Compiling disruption (population affected) per target per storm with {n_proc} threads")
    if n_proc > 1:
        import multiprocessing
        with multiprocessing.Pool(processes=n_proc) as pool:
            disruption_by_target_by_storm: list[tuple[str, pd.DataFrame]] = pool.map(disruption_by_target, snakemake.input.disruption)
        disruption_by_target_by_storm = dict(disruption_by_target_by_storm)
    else:
        disruption_by_target_by_storm = {}
        for disruption_path in snakemake.input.disruption:
            event_id, disruption = disruption_by_target(disruption_path)
            disruption_by_target_by_storm[event_id] = disruption

    # calculate number of years between first and last storm event, necessary for expected annual disruption
    event_ids: list[str] = list(disruption_by_target_by_storm.keys())
    years: set[int] = set(track_year.loc[event_ids, "year"])
    span_years: int = max([1, max(years) - min(years)])  # with a minimum of one

    # target rows, threshold value columns, values are population affected for a given storm
    disruption_all_storms: pd.DataFrame = pd.concat(disruption_by_target_by_storm.values()).groupby("target").sum()

    # create a lookup between target id and the region to which the target's representative point lies within
    logging.info("Creating target to region mapping")
    # filter out targets that are never exposed, we don't need to do an expensive sjoin on them
    target_rep_points: gpd.GeoDataFrame = targets.loc[disruption_all_storms.index, ["population", "geometry"]].copy()
    target_rep_points.geometry = target_rep_points.geometry.representative_point()
    target_to_region_mapping: pd.DataFrame = target_rep_points.sjoin(regions, how="left").drop(columns=["geometry", "index_right"])

    # merge with regions and sum targets across regions
    disruption_by_region = \
        target_to_region_mapping.drop(columns=[f"NAME_{admin_level}"]).merge(disruption_all_storms, on="target").groupby(f"GID_{admin_level}").sum()

    # take the disruption counts and divide by the years passing between first and last storm
    # this division is aligned on the indicies (both set to target ids)
    logging.info("Calculating expected annual population disrupted")
    disruption_fraction_by_region = \
        disruption_by_region.drop(columns=["population"]).divide(span_years, axis=0)

    # merge geometry and name columns back in
    disruption_with_geometry = \
        disruption_fraction_by_region.merge(regions[[f"NAME_{admin_level}", f"GID_{admin_level}", "geometry"]], on=f"GID_{admin_level}", how="right")
    # merge nominal lengths by region back in, too
    disruption_with_population = disruption_with_geometry.merge(disruption_by_region[["population"]], on=f"GID_{admin_level}")

    # write out to disk
    logging.info("Writing out with region geometry")
    gpd.GeoDataFrame(disruption_with_population).to_parquet(snakemake.output.total_disruption_by_region)
