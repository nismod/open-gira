# Extracts

The bounding boxes in calculated from the datasets [earlier](bbox.md)
are chopped up into smaller bounding boxes based on the `slice_count` asked
for in the [config file](../configuration.md).
The `./results/json/tanzania-latest.json` file that defines the bounding box for the
`tanania-latest` dataset is subdivided into a grid of `slice_count` smaller bounding boxes.
The output is saved as `./results/json/tanzania-latest_extracts.geojson`.
