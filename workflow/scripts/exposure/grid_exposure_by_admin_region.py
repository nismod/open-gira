"""
For a given country and storm set, combine all the storm exposure estimates and
aggregate them to some administrative level.

N.B. The reading of exposure files is parallelised with multiprocessing, degree
of parallelism set by calling rule's threads setting.

Write out expected annual exposure (fraction of region's grid edges total
length exposed to winds in excess of a threshold).
"""


import logging
import os
import sys

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr


def plot_event_distributions(thresholds: list[str], exposure_by_event: pd.DataFrame, plot_dir: str) -> None:
    """
    Given exposure data (lengths of edge exposed to wind speed in excess of a
    given threshold), plot the distribution of events.

    Args:
        thresholds: Wind speed thresholds for exposure data, used for indexing
            column from `exposure_by_event`.
        exposure_by_event: Lengths exposed in meters, with event index, threshold columns.
        plot_dir: Location to save histograms.
    """
    # find extrema and number of bins based on all thresholds, so x-axis is
    # constant between thresholds
    x_min_log10 = np.log10(max([1000, 0.1 * min(exposure_by_event.min())]))
    x_max_log10 = np.log10(max(exposure_by_event.max()))
    n_bins = max([10, int(np.ceil(np.cbrt(len(exposure_by_event))))])
    bins = np.logspace(x_min_log10, x_max_log10, n_bins)

    # find y (frequency) maxima across thresholds
    y_max = 0
    for threshold in thresholds:
        frequency, bin_edges = np.histogram(exposure_by_event.loc[:, threshold], bins=bins)
        y_max = max([y_max, max(frequency)])
    y_max *= 1.1

    for threshold in thresholds:

        # filter out zero values
        data = exposure_by_event[threshold].copy()
        non_zero_data = data[data > 0]

        f, ax = plt.subplots()

        ax.hist(
            non_zero_data,
            bins=bins,
            color="green",
            alpha=0.5,
            label="Exposure distribution"
        )
        p90 = np.quantile(non_zero_data, 0.9)
        p95 = np.quantile(non_zero_data, 0.95)
        p99 = np.quantile(non_zero_data, 0.99)
        ax.axvline(p90, label=r"$p_{90}$ = " + f"{p90:.2E}", ls="--", color="pink", alpha=0.7)
        ax.axvline(p95, label=r"$p_{95}$ = " + f"{p95:.2E}", ls="--", color="purple", alpha=0.7)
        ax.axvline(p99, label=r"$p_{99}$ = " + f"{p99:.2E}", ls="--", color="red", alpha=0.7)
        ax.set_ylim(0, y_max)
        ax.set_xlim(10**x_min_log10, 5 * 10**x_max_log10)
        ax.set_xscale("log")
        ax.set_xlabel("Region exposure [m]")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Exposure event distribution, n={len(exposure_by_event)} events\nExposed edge length @ {threshold} [m s-1] threshold")
        ax.grid(alpha=0.1)
        ax.legend()

        f.savefig(os.path.join(plot_dir, threshold.replace(".", "p") + ".png"))

    return


def exposure_by_edge(exposure_path: str) -> pd.DataFrame:
    """
    Given an exposure file, read netCDF and transform to pandas dataframe.

    Args:
        exposure_path: Path to an exposure netCDF containing `length_m`
            variable on `event_id` (singleton), `threshold` and `edge`
            coordinates.

    Returns:
        Exposure dataframe with edge rows and threshold columns for given
            storm. Additional column of `event_id`. Values are lengths (m) of edge
            exposed to wind speeds in excess of threshold speed.
    """

    exposure = xr.open_dataset(exposure_path)

    event_id: str = exposure.event_id.item()
    logging.info(event_id)

    df: pd.DataFrame = exposure.length_m.to_pandas().T

    # N.B. arrow specification requires string column names (but threshold column names are float)
    df.columns = df.columns.astype(str)

    # repeat the event id for each edge
    df["event_id"] = event_id

    return df


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # load network edges
    logging.info("Loading edges")
    edges: gpd.GeoDataFrame = gpd.read_parquet(snakemake.input.grid_edges)
    if edges.empty is True:
        logging.info("No grid representation, write out empty exposure")
        gpd.GeoDataFrame({"geometry": []}, crs=4326).to_parquet(snakemake.output.total_exposure_by_region)
        sys.exit(0)

    edges = edges[["id", "geometry"]]
    edges = edges.rename(columns={"id": "edge"}).set_index("edge")
    edges = edges.set_crs(epsg=4326)
    # calculate length of each edge (transmission and distribution lines)
    edges["nominal_length_m"] = edges.to_crs(edges.estimate_utm_crs()).geometry.length

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

    n_proc: int = snakemake.threads
    logging.info(f"Compiling exposure lengths per edge per storm with {n_proc} threads")
    if n_proc > 1:
        import multiprocessing
        with multiprocessing.Pool(processes=n_proc) as pool:
            exposure_by_edge_by_storm: list[pd.DataFrame] = pool.map(exposure_by_edge, snakemake.input.exposure)
    else:
        exposure_by_edge_by_storm = []
        for exposure_path in snakemake.input.exposure:
            exposure_by_edge_by_storm.append(exposure_by_edge(exposure_path))

    # storm-edge rows (repeated edges), threshold value columns, values are exposed length in meters for a given storm
    exposure_all_storms: pd.DataFrame = pd.concat(exposure_by_edge_by_storm)

    # calculate number of years between first and last storm event, necessary for expected annual exposure
    event_ids: list[str] = list(set(exposure_all_storms.event_id))
    years: set[int] = set(track_year.loc[event_ids, "year"])
    span_years: int = max([1, max(years) - min(years)])  # with a minimum of one

    # create a lookup between edge id and the region to which the edge's representative point lies within
    logging.info("Creating edge to region mapping")
    # filter out edges that are never exposed, we don't need to do an expensive sjoin on them
    edge_rep_points: gpd.GeoDataFrame = edges.loc[edges.reset_index().edge.isin(exposure_all_storms.index)].copy()
    edge_rep_points.geometry = edge_rep_points.geometry.representative_point()
    edge_to_region_mapping: pd.DataFrame = edge_rep_points.sjoin(regions, how="left").drop(columns=["geometry", "index_right"])

    # whole country event distributions
    per_event = exposure_all_storms.groupby("event_id").sum()
    os.makedirs(snakemake.output.country_event_distributions)
    plot_event_distributions(list(per_event.columns), per_event, snakemake.output.country_event_distributions)

    # TODO: per region event distributions?
    #per_event.to_parquet(snakemake.input.per_event_exposure)

    exposure_all_storms = exposure_all_storms.reset_index().drop(columns="event_id").set_index("edge")

    logging.info("Aggregating to region level")
    # edge rows, threshold value columns, values are exposed length in meters as a result of all storms
    exposure_total = exposure_all_storms.groupby("edge").sum()
    # merge with regions and sum edges across regions
    exposure_by_region = \
        edge_to_region_mapping.drop(columns=[f"NAME_{admin_level}"]).merge(exposure_total, on="edge").groupby(f"GID_{admin_level}").sum()
    # take the exposure lengths and divide by the product of original, undamaged lengths and years passing between first and last storm
    # this division is aligned on the indicies (both set to edge ids)
    # we now have an expected annual exposure
    logging.info("Calculating expected annual exposure")
    exposure_fraction_by_region = \
        exposure_by_region.drop(columns=["nominal_length_m"]).divide(exposure_by_region["nominal_length_m"] * span_years, axis=0)
    # merge geometry and name columns back in
    exposure_with_geometry = \
        exposure_fraction_by_region.merge(regions[[f"NAME_{admin_level}", f"GID_{admin_level}", "geometry"]], on=f"GID_{admin_level}", how="right")
    # merge nominal lengths by region back in, too
    exposure_with_length = exposure_with_geometry.merge(exposure_by_region[["nominal_length_m"]], on=f"GID_{admin_level}")
    # write out to disk
    logging.info("Writing out with region geometry")
    gpd.GeoDataFrame(exposure_with_length).to_parquet(snakemake.output.total_exposure_by_region)
