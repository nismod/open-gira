# Rail

The process for creating a rail network is essentially the same as for road.

We can create a topologically connected rail network for a given area from
[OpenStreetMap](https://www.openstreetmap.org) (OSM) data. The resulting
network can be annotated with data retrieved from OSM, along with data looked up
from user-supplied sources (e.g. rehabilitation costs). The network edges will
be labelled with from nodes and to nodes, describing the connectedness of the
network.

## Description

1. The target OSM datasets are downloaded or copied and saved as
   `<output_dir>/input/<dataset>.osm.pbf`.
1. The initial OSM datasets are filtered, keeping only relevant tags
   (using `osmium tags-filter`). This results in smaller files
   `<output_dir>/input/<dataset>_filter-rail.osm.pbf`, where `<dataset>` is the
   key name. 
1. The OSM dataset's headers are examined for a `bbox` property and that is used
   to determine the bounding box for the whole area (`<output_dir>/json/<dataset>.json`).
1. The OSM dataset bounding box is sliced into a grid of smaller bounding boxes
   according to the `slice_count` config option, and these slices are saved
   in a json file `<output_dir>/json/<dataset>-extracts.geojson`.
1. The filtered OSM file is sliced into areas of equal size using the bounding
   box grid. The slices are saved to
   `<output_dir>/slices/<dataset>_filter-rail/slice-<N>.osm.pbf`.
1. Each filtered OSM dataset slice is then converted to the GeoParquet data format,
   resulting in `<output_dir>/geoparquet/<dataset>_filter-rail/raw/slice-<N>_edges.geoparquet`.
1. These geoparquet slices are annotated with various attributes and processed into connected networks
   `<output_dir>/geoparquet/<dataset>_filter-rail/processed/slice-<N>_edges.geoparquet`
1. The slices are joined together into a complete connected network
   `<output_dir>/<dataset>_filter-rail/edges.geoparquet`

## Configuration

To specify a desired network:
- Review and amend the spreadsheets in `bundled_data/transport`, these supply
  information that is used to gap-fill or extend what can be determined from OSM alone.
- Review and amend `config/config.yaml`:
    - The `infrastructure_datasets` map should contain a key pointing to an `.osm.pbf`
      file URL for desired area. There are currently entries for the planet,
      for (some definition of) continents and several countries. We use
      the [geofabrik](http://download.geofabrik.de/) service for continent and
      country-level OSM extracts.
    - Check the OSM filter file pointed to by `network_filters.rail`.
      This file specifies which [elements](https://wiki.openstreetmap.org/wiki/Elements)
      (nodes, ways or relations) to keep (or reject) from the multitude of data
      in an OSM file. See the filter expressions section
      [here](https://docs.osmcode.org/osmium/latest/osmium-tags-filter.html)
      for more information on the syntax of these files.
    - Check and amend `keep_tags.rail`. This list of strings specifies which
      `tags` (attributes) to retain on the filtered elements we extract from
      the `.osm.pbf` file.
    - Review `slice_count`. This controls the degree of parallelism possible.
      With it set to 1, there is no spatial slicing (we create the network in
      a single chunk). To speed network creation for large domains, it can be
      set to a larger square number. The first square number greater than your
      number of available CPUs is a good heuristic.
    - Check and amend the values of `transport.rail`, which provide some
      defaults for OSM data gap-filling.

## Running

And to create the network, by way of example:
```bash
snakemake --cores all -- results/egypt-latest_filter-rail/edges.geoparquet
```

Note that the nodes file, `results/egypt-latest_filter-rail/nodes.geoparquet`
will by default contain the stations and their names as recorded in OSM.