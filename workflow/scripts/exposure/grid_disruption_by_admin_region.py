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

    # event rows, threshold value columns, values are population affected summed across all targets
    disruption_by_event: pd.DataFrame = pd.read_parquet(snakemake.input.disruption_by_event)
    if len(disruption_by_event.index) == 0:
        logging.info("No disruption data, write out empty disruption")
        regions.to_parquet(snakemake.output.expected_annual_disruption)
        sys.exit(0)

    # calculate number of years between first and last storm event, necessary for expected annual disruption
    event_ids: list[str] = list(disruption_by_event.index)
    years: set[int] = set(track_year.loc[event_ids, "year"])
    span_years: int = max([1, max(years) - min(years)])  # with a minimum of one
    logging.info(f"Using {len(event_ids):,d} events, over {span_years:,d} years")

    # target rows, threshold value columns, values are population affected summed across all events
    disruption_by_target: pd.DataFrame = pd.read_parquet(snakemake.input.disruption_by_target)

    # create a lookup between target id and the region to which the target's representative point lies within
    logging.info("Creating target to region mapping")
    # filter out targets that are never exposed, we don't need to do an expensive sjoin on them
    target_rep_points: gpd.GeoDataFrame = targets.loc[:, ["population", "geometry"]].copy()
    target_rep_points.geometry = target_rep_points.geometry.representative_point()
    target_to_region_mapping: pd.DataFrame = target_rep_points \
        .sjoin(regions, how="left") \
        .drop(columns=["geometry", "index_right"])
    # rename index column from "id" to "target" in preparation for merge with `disruption_by_target`
    target_to_region_mapping.index = target_to_region_mapping.index.rename("target")

    # merge with regions and sum targets across regions
    # use how="left" to take every region, exposed or not, to give a complete table
    disruption_by_region = \
        target_to_region_mapping \
        .drop(columns=[f"NAME_{admin_level}"]) \
        .merge(disruption_by_target, on="target", how="left") \
        .groupby(f"GID_{admin_level}") \
        .sum()

    # take the disruption counts and divide by the years passing between first and last storm
    # this division is aligned on the indicies (both set to target ids)
    logging.info("Calculating expected annual population disrupted")
    disruption_fraction_by_region = \
        disruption_by_region.drop(columns=["population"]).divide(span_years, axis=0)

    # merge geometry and name columns back in
    disruption_with_geometry = \
        disruption_fraction_by_region \
        .merge(regions[[f"NAME_{admin_level}", f"GID_{admin_level}", "geometry"]], on=f"GID_{admin_level}")

    # write out to disk
    logging.info("Writing out with region geometry")
    gpd.GeoDataFrame(disruption_with_geometry).to_parquet(snakemake.output.expected_annual_disruption)
