# Composite network generation

Create road or rail networks with varying asset filters (to primary routes in
one region, tertiary in others, etc.).

## Description

For each `infrastructure_dataset` and `network_filter`:
1. Download or copy `.osm.pbf` file to input data directory.
1. Filter OSM file with `osmium tags-filter` for elements matching the filter file (see configuration, below).
1. Define a bounding box for the OSM file, write to disk as JSON.
1. Define a grid partitioning the bounding box into a square number of 'slices' as a series of JSON files, one for each slice.
1. Cut the OSM file into these slices.
1. Convert the sliced OSM files into geoparquet, retaining the `keep_tags` as configured.
1. Clean and annotate features in the geoparquet files (joining additional data such as country, rehabiliation costs, etc.).
1. Join sliced network components together.
1. Join all datasets together.

## Configuration

- Review and amend `config/composite_networks/*.csv`
  - Filter files should include two columns, infrastructure_dataset and
    network_filter.
    For example:
    ```bash
    infrastructure_dataset,network_filter
    thailand-latest,road-secondary
    laos-latest,road-primary
    cambodia-latest,road-primary
    myanmar-latest,road-primary
    ```
    These then map to the infrastructure datasets and network filters
    specified in the `config/config.yaml` file.
- Review and amend `config/config.yaml`:
  - The `composite_network` mapping is from an identifier key to a filter file
    path.
  - The `infrastructure_datasets` map should contain a key pointing to an `.osm.pbf`
    file URL for the desired areas.
    There are currently entries for the planet, for (some definition of)
    continents and several countries. We use the
    [geofabrik](http://download.geofabrik.de/) service for continent and
    country-level OSM extracts.
  - Check the OSM filter file pointed to by `network_filters.road`.
    This file specifies which [elements](https://wiki.openstreetmap.org/wiki/Elements)
    (nodes, ways or relations) to keep (or reject) from the multitude of data
    in an OSM file. See the filter expressions section
    [here](https://docs.osmcode.org/osmium/latest/osmium-tags-filter.html)
    for more information on the syntax of these files.
  - Check and amend `keep_tags.road` and/or `keep_tags.rail`. This list of
    strings specifies which `tags` (attributes) to retain on the filtered
    elements we extract from the `.osm.pbf` file.
  - Review `slice_count`. This controls the degree of parallelism possible.
    With it set to 1, there is no spatial slicing (we create the network in
    a single chunk). To speed network creation for large domains, it can be
    set to a larger square number. The first square number greater than your
    number of available CPUs is a good heuristic.

## Creation

Here's an example creation command:
```bash
snakemake --cores all results/composite_network/south-east-asia-road/edges.gpq
```
