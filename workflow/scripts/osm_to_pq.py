"""OSM.pbf to parquet using osmium, geopandas
"""
import logging
import os
import sys

import geopandas
import osmium
import shapely.wkb as wkblib

# A global factory that creates WKB from a osmium geometry
wkbfab = osmium.geom.WKBFactory()

class WayHandler(osmium.SimpleHandler):
    def __init__(self):
        osmium.SimpleHandler.__init__(self)
        self.output_data = []

    def way(self, w):
        try:
            wkb = wkbfab.create_linestring(w)
            line = wkblib.loads(wkb, hex=True)
        except RuntimeError as e:
            print("Ignoring", e)
            return
        highway = w.tags['highway'] if 'highway' in w.tags else None
        id = w.tags['id'] if 'id' in w.tags else None
        feature = {
            'geometry': line,
            'highway': highway,
            'id': id,
        }
        self.output_data.append(feature)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    logging.info("Start")
    try:
        pbf_path = snakemake.input[0]
        output_path = snakemake.output[0]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        pbf_path, output_path = sys.argv[1:]

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

    print(f"Converting {pbf_path} to .geoparquet")

    h = WayHandler()
    h.apply_file(pbf_path, locations=True)
    geopandas.GeoDataFrame(h.output_data).to_parquet(
        output_path
    )
    # print(geopandas.read_parquet(output_path))
