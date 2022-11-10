# OpenStreetMap data

The OpenStreetMap data is defined in the [config file](../configuration.md) as
a named list of .osm.pbf files.
Each one of these is downloaded (if it's located on the internet) or copied
(if it's located elsewhere) and saved using the name given.

In this case, `https://download.geofabrik.de/africa/tanzania-latest.osm.pbf`
is downloaded and saved as `./results/input/tanzania-latest.osm.pbf`,
and `https://download.geofabrik.de/europe/great-britain/wales-latest.osm.pbf`
is downloaded and saved as `./results/input/wales-latest.osm.pbf`.