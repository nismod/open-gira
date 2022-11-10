# Bounding box calculation

The workflow will slice up the 
`./results/input/<dataset>_filter-highway-core.osm.pbf` 
files into smaller slices to speed up the workflow.
The first step to doing this is to determine the bounding box for the entire .osm.pbf
file, and this is done using `osmium` to read the `fileinfo`.
When this is done, the resulting output is saved as `./results/json/<dataset>.json`.

The bounding box information will be used later in two places.
- When hazard raster files are [trimmed](trim-hazard.md) to match the .osm.pbf files
- When the .osm.pbf files are [sliced](extracts.md) to increase workflow throughput
