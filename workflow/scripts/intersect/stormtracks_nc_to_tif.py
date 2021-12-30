"""Preprocess NetCDF tropical cyclone fixed return period files
- fix metadata to comply with CF conventions
- write TIFFs for relative convenience - matches other hazard layers
"""
import os
from glob import glob

import rioxarray
import xarray
import sys
import json
import ast

# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

    region = "WP"

else:  # linux
    region = sys.argv[1]


def main():
    storm_path = os.path.join("data", "stormtracks", "fixed")
    # fix metadata
    for fname in get_fnames(storm_path):
        with xarray.open_dataset(fname) as ds:
            ds = add_meta(ds)
            for var in list(ds):
                for rp in list(ds["mean"].rp.data):
                    out_fname = fname.replace(".nc", f"_rp{int(rp)}_{var}.tif")
                    out_fname = out_fname.replace("fixed", os.path.join("fixed", "extracted"))
                    print(out_fname)
                    ds[var].sel(rp=rp).rio.to_raster(out_fname)



def get_fnames(storm_path):
    return [
        fname for fname in glob(os.path.join(storm_path, "*.nc")) if f"RETURN_PERIODS_{region}" in fname and "withmeta" not in fname
    ]


def add_meta(ds):
    ds = ds.transpose("rp", "lat", "lon")
    ds.coords["lon"].attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "axis": "X",
    }
    ds.rio.write_crs("epsg:4326", inplace=True)
    return ds

main()
#if __name__ == "__main__":
#    main()