# Split world

The [input data](../download/power_download.md) consists of large files which is not possible to process as one for the purposes of this analysis. The solution is to therefore split the world into boxes (each with their own unique box_id) and, where possible, process each box individually. This has advantages of parallel processing with the snakemake workflow too.

### Specifying the box dimensions
Each box is a square with width and height equal to the `box_width_height` value (units in degrees) specified in the config.yaml file. The result of this is depicted by the following image (`box_width_height=5`).
![Example world split `box_width_height=5`](../power_img/world_splitter.png)

### Specifying individual boxes
Futhermore, the option to process and intersect within a set of boxes is possible. In the config.yaml file, a list of boxes may be specified for the `specific_boxes` parameter. This means the analysis is now not global and there may be storms which do not intersect with these boxes. Additionally, it must be known which box box_id values are required in advance which can be done through geometrical equations or, more simply, by running split world first and then visualising `results/power_process/world_boxes.gpkg` and `results/input/adminboundaries/gadm36_levels.gpkg` (level 0 recommended) in QGIS (for example). 