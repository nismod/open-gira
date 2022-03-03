# Intersect infrastructure and flood level

The intersection process inspects the geographic data for each of the rows in our .geoparquet files in 
`./results/geoparquet/<dataset>_filter-<filter>/` and produces a new file with the flood depths affecting that road in each
of the flood scenarios in `./results/input/hazard-<hazard>/<dataset>/`.
These new files are saved in `./results/splits/<dataset>_filter-<filter>/hazard-<hazard>/`. 
Let's use the Python console or update our Python script to load in the intersection information for
slice 32 of Tanzania:

```python
# We already imported geopandas
slice_file_name = 'results/splits/tanzania-latest_filter-highway-core/hazard-aqueduct-river/slice-32.geoparquet'
slice_gp = geopandas.read_parquet(slice_file_name)
print(slice_gp)
```

We should see something like:
```text
                                              geometry  ...  inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100
0    LINESTRING (39.53152 -6.26680, 39.53157 -6.266...  ...                                           0.0
0    LINESTRING (39.53176 -6.26712, 39.53185 -6.267...  ...                                           0.0
0    LINESTRING (39.53359 -6.27546, 39.53375 -6.275...  ...                                           0.0
0    LINESTRING (39.53594 -6.28379, 39.53595 -6.283...  ...                                           0.0
0    LINESTRING (39.53644 -6.28704, 39.53660 -6.288...  ...                                       -9999.0
..                                                 ...  ...                                           ...
855  LINESTRING (39.28642 -6.79626, 39.28654 -6.796...  ...                                       -9999.0
856  LINESTRING (39.28679 -6.79718, 39.28696 -6.797...  ...                                       -9999.0
856  LINESTRING (39.28847 -6.80050, 39.28849 -6.800...  ...                                       -9999.0
857  LINESTRING (39.28782 -6.80252, 39.28790 -6.802...  ...                                       -9999.0
857  LINESTRING (39.28841 -6.80050, 39.28838 -6.800...  ...                                       -9999.0
[2504 rows x 15 columns]
```

Notice that there are more rows and columns than there were before.
There are more rows because each highway segment has been split up into one row for each cell in the raster grid
that it passes through.[^raster]
There are more columns because each hazard scenario we had in 
`./results/input/hazard-<hazard>/<dataset>/` has added a column 
with its filename as the colum name and the maximum flood depth for 
each stretch of highway in each raster cell as its values.

Corresponding `.parquet` files (without geometries) are also created.

[^raster]: The flood maps in `./results/input/hazard-<hazard>/<dataset>/` are raster files that show flood depth in each cell.