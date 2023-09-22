"""
Concatenate exposure files (per-storm) into a single file.
"""

import logging

import pandas as pd
import numpy as np
import xarray as xr


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

    n_proc: int = snakemake.threads
    logging.info(f"Compiling exposure lengths per edge per storm with {n_proc} threads")
    if n_proc > 1:
        import multiprocessing
        with multiprocessing.Pool(processes=n_proc) as pool:
            exposure_by_edge_by_storm: list[pd.DataFrame] = pool.map(exposure_by_edge, snakemake.input.exposure_by_event)
    else:
        exposure_by_edge_by_storm = []
        for exposure_path in snakemake.input.exposure_by_event:
            exposure_by_edge_by_storm.append(exposure_by_edge(exposure_path))

    # storm-edge rows (repeated edges), threshold value columns, values are exposed length in meters for a given storm
    logging.info("Concatenating exposure by storm")
    # this step is very memory intensive, ~ 100GB for STORM-constant and the Philippines
    exposure_all_storms: pd.DataFrame = pd.concat(exposure_by_edge_by_storm)

    logging.info("Writing to disk")
    exposure_all_storms.to_parquet(snakemake.output.concatenated)
