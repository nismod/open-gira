# Filter each slice by infrastructure

The original .osm.pbf files (`./results/input/<dataset>.osm.pbf`) 
contain lots and lots of information about points of interest, waterways, etc.
We strip them down to versions that contain only information about highways 
(and only highways that match one of the types we specified in 
`./config/highway-core.txt`, or whichever other file we declared
as the `osmium_tags_filters_file` variable in `./config/config.yaml`).

This file is in `./results/input`, so let's open it up in QGIS.
It will have been renamed to indicate that it's been filtered 
(using the filename of the `osmium_tags_filters_file`).
Load `./results/input/tanzania-latest_filter-highway-core.osm.pbf`, 
and notice that there are fewer roads than in the previous file.
We've coloured the roads in the new file in red.

![QGIS screenshot showing red roads overlaid on black ones.](../../img/QGIS-filtered.png)

If you zoom in a lot, you'll see there are many, many roads that are not included in this smaller file.
Here, we've zoomed in on the port area of Dar es Salaam.

![QGIS screenshot showing many black roads, a few of which have red roads on top.](../../img/QGIS-filtered_zoom.png)

If you want to, you can do the same for the Wales datasets.