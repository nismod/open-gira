# Convert Osmium files to GeoParquet format

Geoparquet files are files that contain geographic data alongside other information.
They are loaded into Python as dataframes (`geopandas.GeoDataFrame`).
The .geoparquet files are stored in `./results/geoparquet`.
Once again, let's open up the slice for Dar es Salaam, slice 32.
We'll do this using Python, so either write a script or pop open the console, and put in the following:

```python
import geopandas  # installed when we installed open-gira

# Assuming the open-gira root is the current working directory:
file_name = 'results/geoparquet/tanzania-latest_filter-highway-core_slice-32.geoparquet'

gp = geopandas.read_parquet(file_name)  # load the file

print(gp)  # Take a peek at the dataframe
```

We should get an output something like:

```text
                                              geometry     highway    id
0    LINESTRING (39.13700 -6.62114, 39.13776 -6.621...       trunk  None
1    LINESTRING (39.26407 -6.77160, 39.26266 -6.771...   secondary  None
2    LINESTRING (39.23765 -6.77090, 39.23705 -6.771...   secondary  None
3    LINESTRING (39.30578 -6.82119, 39.30526 -6.821...   secondary  None
4    LINESTRING (39.37496 -6.86819, 39.37437 -6.868...   secondary  None
..                                                 ...         ...   ...
865  LINESTRING (39.24408 -6.84490, 39.24324 -6.845...       trunk  None
866  LINESTRING (39.24554 -6.84430, 39.24436 -6.844...  trunk_link  None
867  LINESTRING (39.24168 -6.84576, 39.24211 -6.845...  trunk_link  None
868  LINESTRING (39.24986 -6.84211, 39.24903 -6.842...  trunk_link  None
869  LINESTRING (39.24591 -6.84376, 39.24603 -6.843...  trunk_link  None

[870 rows x 3 columns]
```

That looks about right -- we can see we have some geometry information, the type of the highway the geometry
describes, and an (unused) id field.
