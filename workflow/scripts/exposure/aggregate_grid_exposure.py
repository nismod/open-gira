"""
Take exposure files (one per storm, with per-edge data) and aggregate into sample per-event total and sample per-edge total.
"""

import glob
import logging
import os

import dask
import dask.dataframe
import pandas as pd
import xarray as xr


def make_schema(column_dtypes: dict[str, type]) -> pd.DataFrame:
    """
    Make an empty dataframe with column names and dtypes.

    Args:
        column_dtypes: Mapping from column name to type of column.

    Returns:
        Empty dataframe with column names and data types set.
    """
    schema = pd.DataFrame(index=None)
    for name, dtype in column_dtypes.items():
        schema[name] = pd.Series(dtype=dtype)
    return schema


@dask.delayed
def exposure_by_edge(exposure_path: str, schema: pd.DataFrame) -> pd.DataFrame:
    """
    Given an exposure file, read netCDF and transform to pandas dataframe.

    Args:
        exposure_path: Path to an exposure netCDF containing `length_m`
            variable on `event_id` (singleton), `threshold` and `edge`
            coordinates.
        schema: Empty dataframe with schema (column names and dtypes) of
            desired output

    Returns:
        Exposure dataframe with edge rows and threshold columns for given
            storm. Additional column of `event_id`. Values are lengths (m) of edge
            exposed to wind speeds in excess of threshold speed.
    """

    exposure = xr.open_dataset(exposure_path)

    event_id: str = exposure.event_id.item()

    df: pd.DataFrame = exposure.length_m.to_pandas().T

    # N.B. arrow specification requires string column names (but threshold column names are float)
    df.columns = [f"{value:.1f}" for value in df.columns.values]

    # take the netCDF edge dimension (currently the dataframe index) and store as data
    df: pd.DataFrame = df.reset_index()

    # repeat the event id for each edge
    df["event_id"] = event_id

    # in the case of missing columns, add them by concatenating with schema
    df: pd.DataFrame = pd.concat([schema, df])

    return df


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Compiling exposure lengths per edge per storm")

    # we're using dask so we can concatenate on disk rather than in memory
    # (pd.concat is a disaster in terms of memory usage)
    # for 10,000 years of storms impacting Japan, a single thread will take ~7min
    dask.config.set(scheduler='single-threaded')

    # construct a dataframe with no data, but the appropriate columns and dtypes
    # exposure_by_edge will coerce the data to match this schema
    threshold_cols: list[str] = [f"{value:.1f}" for value in snakemake.params.thresholds]
    column_dtypes: dict[str, type] = {"event_id": str, "edge": int, **dict(zip(threshold_cols, [float] * len(threshold_cols)))}
    schema: pd.DataFrame = make_schema(column_dtypes)

    # create list of data read and transform tasks
    tasks = []
    for file_path in sorted(glob.glob(f"{snakemake.input.exposure_by_event}/*.nc")):
        tasks.append(exposure_by_edge(file_path, schema))

    # write out to disk as parquet
    if tasks:
        exposure_all_storms = dask.dataframe.from_delayed(tasks)
        exposure_all_storms.drop(columns=["edge"]).groupby("event_id").sum().to_parquet(snakemake.output.by_event)
        exposure_all_storms.drop(columns=["event_id"]).groupby("edge").sum().to_parquet(snakemake.output.by_edge)
    else:
        # we have declared to snakemake that this output will be a directory (when there's data, it's sharded)
        # write out the schema as a pretend directory parquet dataset
        os.makedirs(snakemake.output.by_event)
        os.makedirs(snakemake.output.by_edge)
        schema.drop(columns=["edge"]).groupby("event_id").sum().to_parquet(os.path.join(snakemake.output.by_event, "part.0.parquet"))
        schema.drop(columns=["event_id"]).groupby("edge").sum().to_parquet(os.path.join(snakemake.output.by_edge, "part.0.parquet"))