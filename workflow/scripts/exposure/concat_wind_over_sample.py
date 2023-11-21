import logging
import sys

import dask
from dask.delayed import Delayed
import xarray as xr

from open_gira.io import netcdf_packing_parameters
from open_gira.wind import empty_wind_da


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Reading wind fields from each sample")
    # concatenated xarray dataset is chunked by input file
    # N.B. this is lazily loaded
    # cannot use open_mfdataset here, as event_id coordinates are not monotonically increasing
    all_samples = xr.concat(
        [
            xr.open_dataset(path, chunks={"max_wind_speed": 1}).sortby("event_id")
            for path in snakemake.input.sample_paths
        ],
        dim="event_id"
    ).sortby("event_id")

    if all_samples.event_id.size == 0:
        logging.info("Input data empty, writing empty file to disk")
        # write empty netcdf (with appropriate schema) and exit
        empty_wind_da().to_netcdf(snakemake.output.concat)
        sys.exit(0)

    logging.info("Computing packing factors for all samples")
    # we used dask to allow for chunked calculation of data min/max and to stream to disk
    # use the synchronous scheduler to limit dask to a single process (reducing memory footprint)
    scheduler = "synchronous"

    # compute packing factors for all data, need global min and max
    # the implementation below reads all the data chunks twice, once for min and once for max
    # would be nice to avoid this duplication
    scale_factor, add_offset, fill_value = netcdf_packing_parameters(
        all_samples.max_wind_speed.min().compute(scheduler=scheduler).item(),
        all_samples.max_wind_speed.max().compute(scheduler=scheduler).item(),
        16
    )

    logging.info("Writing pooled wind fields to disk")
    serialisation_task: Delayed = all_samples.to_netcdf(
        snakemake.output.concat,
        encoding={
            "max_wind_speed": {
                'dtype': 'int16',
                'scale_factor': scale_factor,
                'add_offset': add_offset,
                '_FillValue': fill_value
            }
        },
        compute=False
    )
    serialisation_task.compute(scheduler=scheduler)
