osmium tags-filter osm/africa-latest.osm.pbf w/highway=motorway,motorway_link,trunk,trunk_link,primary,primary_link,secondary,secondary_link,tertiary,tertiary_link -o osm/africa-latest-highway-core.osm.pbf
python osm_to_pq.py osm/africa-latest-highway-core.osm.pbf osm/
python network_hazard_intersection.py osm/africa-latest-highway-core.geoparquet aqueduct/ outputs/
