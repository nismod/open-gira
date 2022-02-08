# Bounding box calculation

The workflow will slice up the `./data/tanzania-latest_filter-highway-core.osm.pbf` file into
smaller slices to speed up the workflow.
The first step to doing this is to determine the bounding box for the entire .osm.pbf
file, and this is done using `osmium` to read the `fileinfo`.
When this is done, the resulting output is saved as `./results/json/tanzania-latest.json`
and passed onto the step that calculates the slices.