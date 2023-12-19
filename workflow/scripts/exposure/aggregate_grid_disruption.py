"""
Take per-event disruption files with per-target rows and aggregate into a per-target file and a per-event file.
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
def disruption_by_target(file_path: str, schema: pd.DataFrame) -> pd.DataFrame:
    """
    Given a disruption file, read netCDF and transform to pandas dataframe.

    Args:
        file_path: Path to a disruption netCDF containing `customers_affected`
            variable on `event_id` (singleton), `threshold` and `target`
            coordinates.
        schema: Empty dataframe with schema (column names and dtypes) of
            desired output

    Returns:
        Disruption dataframe with target rows and threshold columns for given
            storm. Additional column of `event_id`. Values are numbers of people
            assessed to be at risk of disruption in a given target.
    """

    ds = xr.open_dataset(file_path)

    event_id: str = ds.event_id.item()

    df: pd.DataFrame = ds.customers_affected.to_pandas().T

    # N.B. arrow specification requires string column names (but threshold column names are float)
    df.columns = [f"{value:.1f}" for value in df.columns.values]

    # take the netCDF target dimension (currently the dataframe index) and store as data
    df: pd.DataFrame = df.reset_index()

    # repeat the event id for each target
    df["event_id"] = event_id

    # in the case of missing columns, add them by concatenating with schema
    df: pd.DataFrame = pd.concat([schema, df])

    return df


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Compiling disruption per target per storm")

    # we're using dask so we can concatenate on disk rather than in memory
    # (pd.concat is a disaster in terms of memory usage)
    # for 10,000 years of storms impacting Japan, a single thread will take ~7min
    dask.config.set(scheduler='single-threaded')

    # construct a dataframe with no data, but the appropriate columns and dtypes for output
    # disruption_by_target will coerce the data to match this schema
    threshold_cols: list[str] = [f"{value:.1f}" for value in snakemake.params.thresholds]
    column_dtypes: dict[str, type] = {"event_id": str, "target": int, **dict(zip(threshold_cols, [float] * len(threshold_cols)))}
    schema: pd.DataFrame = make_schema(column_dtypes)

    # create list of data read and transform tasks
    tasks = []
    for file_path in sorted(glob.glob(f"{snakemake.input.disruption_by_event}/*.nc")):
        tasks.append(disruption_by_target(file_path, schema))

    # write out to disk as parquet
    if tasks:
        disruption_all_storms = dask.dataframe.from_delayed(tasks)
        disruption_all_storms.drop(columns=["target"]).groupby("event_id").sum().to_parquet(snakemake.output.by_event)
        disruption_all_storms.drop(columns=["event_id"]).groupby("target").sum().to_parquet(snakemake.output.by_target)
    else:
        # we have declared to snakemake that this output will be a directory (when there's data, it's sharded)
        # write out the schema as a pretend directory parquet dataset
        os.makedirs(snakemake.output.by_event)
        os.makedirs(snakemake.output.by_target)
        schema.drop(columns=["target"]).groupby("event_id").sum().to_parquet(os.path.join(snakemake.output.by_event, "part.0.parquet"))
        schema.drop(columns=["event_id"]).groupby("target").sum().to_parquet(os.path.join(snakemake.output.by_target, "part.0.parquet"))