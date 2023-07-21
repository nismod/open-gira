# Transport - flooding

The pipeline consists in the following steps:

1. The target OSM datasets are downloaded or copied and saved as
   `<output_dir>/input/<dataset>.osm.pbf`.
2. The initial OSM datasets are filtered, keeping only relevant tags for road links
   (using `osmium tags-filter`). This results in smaller files
   `<output_dir>/input/<dataset>_filter-<filters>.osm.pbf`, where `<dataset>` is the
   key name and `<filters>` is the filename of the `osmium_tags_filter` file in the config.
3. The OSM dataset's headers are examined for a `bbox` property and that is used
   to determine the bounding box for the whole area (`<output_dir>/json/<dataset>.json`).
4. The hazard raster files for each hazard datasets are located by reading a list
   of their locations from the config. Each of these locations is visited and the
   .tif file downloaded or copied to `<output_dir>/input/hazard-<hazard>/raw/<filename>`
   where `<hazard>` is the keyname in the config and `<filename>` is the file's
   base name.
5. Each hazard raster file is clipped to contain just the hazard data for each dataset.
   These files are stored in `<output_dir>/input/hazard-<hazard>/<dataset>/<filename>`
   where `<dataset>` is the OSM dataset whose bounding box is used for clipping.
6. The OSM dataset bounding box is sliced into a grid of smaller bounding boxes
   according to the `slice_count` config option, and these slices are saved
   in a json file `<output_dir>/json/<dataset>-extracts.geojson`.
7. The filtered OSM file is sliced into areas of equal size using the bounding
   box grid from step 6. The slices are saved to
   `<output_dir>/slices/<dataset>_filter-<filter>/slice-<N>.osm.pbf`.
8. Each filtered OSM dataset slice is then converted to the GeoParquet data format,
   resulting in `<output_dir>/geoparquet/<dataset>_filter-<filters>_slice-<N>.geoparquet`.
9. Each geoparquet slice is intersected against flood level data from the
   hazard datasets. The hazard datasets consist of a collection of
   raster data files. The network/hazard intersection results in data
   `<output_dir>/splits/<dataset>_filter-<filters>_slice-<N>_hazard-<hazard>.geoparquet`
   describing roads split according to the raster grid and associated flood level values.
   A corresponding `parquet` files (without geometries) is also created.
10. Split data is then joined into a unique dataset describing
    infrastructure and associated hazard level values for each combination of
    OSM dataset and hazard dataset. This results in
    `<output_dir>/<dataset>_filter-<filters>_hazard-<hazard>.geoparquet`.

These steps, along with the output produced at each stage,
are described in the subsections of this chapter.

These steps are summarised in the digital acyclic graph for `slice_count: 1`, for just the
`tanzania-latest` infrastructure and `aqueduct-coast` hazard data.

[![DAG of the workflow for
the Tanzania dataset and coast flooding data](./img/DAG-simple.png)](./img/DAG-simple.png)
