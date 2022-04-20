# Slice OSM dataset

The filtered dataset (e.g. `./results/input/tanzania-latest_filter-highway-core.osm.pbf`) 
is sliced into areas of equal size according to the slices defined in 
`./results/json/tanzania-latest_extracts.geojson`.
We'll find these in the `results/slices/tanzania-latest_filter-highway-core` directory.
As with the main .osm.pbf file, we can view this in QGIS. 

Open QGIS and load in the _lines_ for the filtered datafile, 
`./results/input/tanzania-latest_filter-highway-core.osm.pbf`.
When that's loaded, pull up one of the slices we created.
Several of the slices will not have road information, so we won't get any feedback 
in QGIS that we've got any information in there at all. 
So we'll pick one that _does_ have information, slice 32 (Dar es Salaam):
`./results/slices/tanzania-latest_filter-highway-core/slice-32.osm.pbf`.
Load up the _lines_ in that and you should see that it doubles up the road lines 
around Dar es Salaam (Eastern coast of Tanzania).
In the image below, we have recoloured the lines so that the lines for slice 32 are in yellow.

![QGIS screenshot showing roads in red with a subsection overlaid in yellow.](transport_img/QGIS-slice32.png)
