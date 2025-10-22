import sys
from pathlib import Path
from osgeo import gdal


gdal.UseExceptions()
DIRECTORY = sys.argv[1]
for fname in list(Path(DIRECTORY).glob("*tif")):
    # grid = snail.intersection.GridDefinition.from_raster(fname)
    # print(fname.stem, "snail transform is", grid.transform)

    with gdal.Open(str(fname)) as ds:
        t = ds.GetGeoTransform()
        print(fname.stem, "gdal transform is", t)
