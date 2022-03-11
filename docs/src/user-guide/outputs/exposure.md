# Exposure rasters

The .geoparquet file containing all the line and flood data can be used to
generate spatial summary information about the hazard-infrastructure intersection.
In this step, a raster file is produced that shows areas exposed to heavy flooding
graded by the amount of road in them.

The steps involved are:
- filter the .geoparquet to remove all rows with flood depth less than `exposure_threshold` m
  (`exposure_threshold` is defined in the [config file](../configuration.md))
- group road segments by cell index and calculate their total length (in km)
- load the original hazard raster file, e.g.
  `./results/input/hazard-aqueduct-river/tanzania-mini/inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif`
- create a copy of that file at
  `./results/exposure/tanzania-mini_filter-highway-core/hazard-aqueduct-river/exposure_inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif`
- replace the flood info band with a road length info band 
  (0 for any cell with less than `exposure_depth` of flooding)
