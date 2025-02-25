"""Set geoTransform values for Deltares coastal flooding extracts

- NASADEM and MERITDEM tiffs have practically the same but unequal geotransforms
  so override to specify exact floating point values
"""

from pathlib import Path

import snail.intersection
from osgeo import gdal

gdal.UseExceptions()

# 90m resolution values:
# SET_TRANSFORM = (
#     -180.00041666666667,
#     0.0008333333333333332,
#     0.0,
#     72.00041666666667,
#     0.0,
#     -0.0008333333333333334,
# )

# 1km resolution values:
SET_TRANSFORM = (-180.00041666666667, 0.01, 0.0, 72.00958333333332, 0.0, -0.01)

DIRECTORY = "results/input/hazard-coastal-deltares/planet-latest"

paths = list(Path(DIRECTORY).glob("*tif"))
for p in paths:
    with gdal.Open(str(p)) as ds:
        old = ds.GetGeoTransform()
        print(p.stem, "was", old)
        ds.SetGeoTransform(SET_TRANSFORM)


for fname in list(Path(DIRECTORY).glob("*tif")):
    grid = snail.intersection.GridDefinition.from_raster(fname)
    print(fname.stem, "is", grid.transform)
