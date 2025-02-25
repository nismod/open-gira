import sys
from pathlib import Path

from osgeo import gdal

gdal.UseExceptions()

DIRECTORY = sys.argv[1]
SET_TRANSFORM = tuple(float(t) for t in sys.argv[2:8])

print("setting", SET_TRANSFORM)


paths = list(Path(DIRECTORY).glob("*tif"))
for p in paths:
    with gdal.Open(str(p)) as ds:
        old = ds.GetGeoTransform()
        print(p.stem, "was", old)
        ds.SetGeoTransform(SET_TRANSFORM)
