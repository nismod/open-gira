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
    pbf_path, outputs_path = sys.argv[1:]
    slug = os.path.basename(pbf_path).replace(".osm.pbf", "")

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

    h = WayHandler()
    h.apply_file(pbf_path, locations=True)
    geopandas.GeoDataFrame(h.output_data).to_parquet(
        os.path.join(
            outputs_path,
            f'{slug}.geoparquet'))
