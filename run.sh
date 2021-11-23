set -e
set -x

continent=$1

osmium tags-filter ../openstreetmap/${continent}-latest.osm.pbf w/highway=motorway,motorway_link,trunk,trunk_link,primary,primary_link,secondary,secondary_link,tertiary,tertiary_link -o data/osm/${continent}-latest-highway-core.osm.pbf
python osm_to_pq.py data/osm/${continent}-latest-highway-core.osm.pbf data/osm/
python network_hazard_intersection.py data/osm/${continent}-latest-highway-core.geoparquet id data/aqueduct/ data/aqueduct/aqueduct_coastal.csv outputs/
python network_hazard_intersection.py data/osm/${continent}-latest-highway-core.geoparquet id data/aqueduct/ data/aqueduct/aqueduct_river.csv outputs/
