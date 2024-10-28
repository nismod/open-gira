# Rail

The process for creating a rail network is essentially the same as for road.

We can create a topologically connected rail network for a given area from
[OpenStreetMap](https://www.openstreetmap.org) (OSM) data. The resulting
network can be annotated with data retrieved from OSM, along with data looked up
from user-supplied sources (e.g. rehabilitation costs).

## Description

1. Download or copy `.osm.pbf` file to input data directory.
1. Filter OSM file with `osmium tags-filter` for elements matching the filter file (see configuration, below).
1. Define a bounding box for the OSM file, write to disk as JSON.
1. Define a grid partitioning the bounding box into a square number of 'slices' as a series of JSON files, one for each slice.
1. Cut the OSM file into these slices.
1. Convert the sliced OSM files into geoparquet, retaining the `keep_tags` as configured.
1. Clean and annotate features in the geoparquet files (joining additional data such as country, rehabiliation costs, etc.).
1. Join sliced network components together.

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

## Creation

And to create the network, by way of example:

```bash
snakemake --cores all -- results/egypt-latest_filter-rail/edges.gpq
```

Note that the nodes file, `results/egypt-latest_filter-rail/nodes.gpq`
will by default contain the stations and their names as recorded in OSM.
