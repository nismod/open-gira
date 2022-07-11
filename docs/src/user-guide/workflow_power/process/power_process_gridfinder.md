# Process gridfinder

This file takes the [gridfinder input data](../download/power_download_gridfinder.md) and assigns each component
(network edge) a box_id, based on its location. Here, the `idxbox()` function from
`workflow/scripts/process/process_power_functions.py` is used to index the components. Through knowledge of the box
parameters from the [world split](power_process_worldsplit.md), it is possible to calculate the box_id:

    box_id = (-np.round(lat / box_width_height + 1 / 2) + lat_max / box_width_height) * num_cols + np.round(lon / box_width_height - 1 / 2) - lon_min / box_width_height

- lat - latitude of component centroid
- lon - longitude of component centroid
- box_width_height - height and width of box
- lat_max - maximum latitutde value
- num_cols - number of cols (i.e. number of boxes across the equator)
- lon_min - minimum longitude value (lon_min can be negative)
- tot_boxes - total number of boxes

The box numbering starts at the top left (north-west), continues left to right (west to east) and
ends at the bottom right (south-east) for longitude and latitude sign convention consistency. The edge cases are
treated separately. Numpy allows for array inputs for lat and lon which improves computational efficiency. Note that if
no components exist in the box_id, then an empty .gpkg file is produced. This is to allow for consistent snakemake
workflows as it can not be known in advance.

