# Convert Osmium files to GeoParquet format

Geoparquet files are files that contain geographic data alongside other information.
They are loaded into Python as dataframes (`geopandas.GeoDataFrame`).
The .geoparquet files are stored in `./results/geoparquet/`.
Once again, let's open up the slice for Dar es Salaam, slice 32.
We'll do this using Python, so either write a script or pop open the console, 
and put in the following:

```python
import geopandas  # installed when we installed open-gira

# Assuming the open-gira root is the current working directory:
file_name = 'results/geoparquet/tanzania-latest_filter-highway-core/slice-32.geoparquet'

gp = geopandas.read_parquet(file_name)  # load the file

print(gp)  # Take a peek at the dataframe
```

We should get an output something like:

```text
                                              geometry  ...  end_node_degree
0    LINESTRING (39.53152 -6.26680, 39.53157 -6.266...  ...                3
1    LINESTRING (39.13700 -6.62114, 39.13776 -6.621...  ...                2
2    LINESTRING (39.26407 -6.77160, 39.26266 -6.771...  ...                2
3    LINESTRING (39.23765 -6.77090, 39.23705 -6.771...  ...                2
4    LINESTRING (39.23340 -6.77427, 39.23306 -6.774...  ...                2
..                                                 ...  ...              ...
853  LINESTRING (39.28411 -6.78375, 39.28413 -6.783...  ...                2
854  LINESTRING (39.28412 -6.78825, 39.28405 -6.787...  ...                2
855  LINESTRING (39.28419 -6.78823, 39.28512 -6.791...  ...                2
856  LINESTRING (39.28679 -6.79718, 39.28696 -6.797...  ...                2
857  LINESTRING (39.28782 -6.80252, 39.28790 -6.802...  ...                2

[858 rows x 12 columns]
```

That looks about right -- we can see we have some geometry information, 
and a bunch of other information. Let's take a look at that other information:

```python
# print out the first 6 rows, without the geometry column
print(gp.loc[:6, gp.columns != 'geometry'].to_string())
```

We should see something like:

```text
     way_id  segment_id tag_highway  start_node_reference  start_node_longitude  start_node_latitude  start_node_degree  end_node_reference  end_node_longitude  end_node_latitude  end_node_degree
0   8412466           0     primary             589993167             39.531516            -6.266800                  1          6120475842           39.564118          -6.410254                3
1  17057163           0       trunk             252499777             39.136997            -6.621140                  2          1409426835           39.140470          -6.622979                2
2  23319781           0   secondary             252498381             39.264073            -6.771600                  2           252498449           39.256254          -6.772167                2
3  23320290           0   secondary             252499782             39.237650            -6.770895                  3           344043335           39.233404          -6.774266                2
4  23320290           1   secondary             344043335             39.233404            -6.774266                  2          1411287465           39.222948          -6.786472                2
5  23321336           0   secondary            3698315617             39.305784            -6.821191                  2          3110649330           39.303854          -6.821512                3
6  23321341           0   secondary             330140099             39.374962            -6.868194                  2          3747757047           39.361199          -6.863670                2
```

There are a bunch of columns there that hold information from the original .osm.pbf file.
Every road in that file has been split up into segments, and each segment contains only a 
start node and an end node that may have junctions.
In other words, we have removed all the middle nodes that do not form junctions, meaning
that we now have a network graph of our roadways.

The columns include the id of the original way in the .osm.pbf file, the segment id of that
way, the tag information that we requested in the `./config/config.yaml` `keep_tags` field,
and then the reference (id), location, and degree of the start and end nodes.
