# Intersect infrastructure and flood level

The intersection process inspects the geographic data for each of the rows in our .geoparquet files in 
`./results/geoparquet/` and produces a new file with the flood depths affecting that road in each
of the flood scenarios in `./data/aqueduct/`.
These new files are saved in `./results/splits/`. 
Let's use the Python console or update our Python script to load in the intersection information for
slice 32:

```python
# We already imported geopandas
slice_file_name = 'results/splits/tanzania-latest_filter-highway-core_slice-32_hazard-aqueduct-river.geoparquet'
slice_gp = geopandas.read_parquet(slice_file_name)
print(slice_gp)
```

We should see something like:
```text
                                               geometry  ... hazard_river__climate_rcp8p5__model_MIROC-ESM-CHEM__y_2080__rp_1000
0     LINESTRING (39.13700 -6.62114, 39.13701 -6.62114)  ...                                                0.0                 
1     LINESTRING (39.13701 -6.62114, 39.13776 -6.621...  ...                                                0.0                 
2     LINESTRING (39.26407 -6.77160, 39.26266 -6.771...  ...                                            -9999.0                 
3     LINESTRING (39.25737 -6.77304, 39.25733 -6.773...  ...                                                0.0                 
4     LINESTRING (39.23765 -6.77090, 39.23740 -6.77120)  ...                                                0.0                 
                                                 ...  ...                                                ...                 
2498  LINESTRING (39.24168 -6.84576, 39.24211 -6.845...  ...                                                0.0                 
2499  LINESTRING (39.24986 -6.84211, 39.24903 -6.842...  ...                                                0.0                 
2500  LINESTRING (39.24903 -6.84267, 39.24808 -6.843...  ...                                                0.0                 
2501  LINESTRING (39.24591 -6.84376, 39.24603 -6.843...  ...                                                0.0                 
2502  LINESTRING (39.24903 -6.84232, 39.24967 -6.84201)  ...                                                0.0     
            
[2503 rows x 283 columns]
```

Notice that there are more rows and columns than there were before.
There are more rows because each highway has been split up into one row for each cell in the raster grid
that it passes through.[^raster]
There are more columns because each hazard scenario we had in `./data/aqueduct` has added a column 
with its filename as the colum name and the maximum flood depth for each stretch of highway in each raster cell
as its values.

Corresponding `.parquet` files (without geometries) are also created.

[^raster]: The flood maps in `./data/aqueduct` are raster files that show flood depth in each cell.