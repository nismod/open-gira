"""
Rasterize the merged protection layer from the FLOPROS shapefile
Going to use the RWI dataset as a reference grid 
"""

import logging

import rasterio
from rasterio.features import rasterize
import fiona

if __name__ == "__main__":

    try:
        flopros_path: str = snakemake.input["flopros"]
        rwi_path: str = snakemake.input["rwi_file"]
        output_path: str = snakemake.output["rasterized_flopros"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Rasterizing the FLOPROS dataset.")

logging.info("Reading reference RWI data.")
with rasterio.open(rwi_path) as ref:
    transform = ref.transform
    width = ref.width
    height = ref.height
    crs = ref.crs
    out_shape = (height, width)

logging.info("Reading the FLOPROS shapefile.")
with fiona.open(flopros_path) as shp:
    shapes = list(
        (feature["geometry"], feature["properties"]["MerL_Riv"])
        for feature in shp
    )

logging.info("Rasterizing FLOPROS.")
raster = rasterize(
    shapes,
    out_shape=out_shape,
    transform=transform,
    fill = 0,
    dtype='float32'
)

logging.info("Writing output raster.")
with rasterio.open(
    output_path,
    "w",
    height=height,
    width=width,
    count=1,
    dtype=raster.dtype,
    crs=crs,
    transform=transform,
) as dst:
    dst.write(raster, 1)

logging.info("Done.")