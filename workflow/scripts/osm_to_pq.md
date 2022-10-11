```
>>> df = gpd.read_parquet("results/geoparquet/egypt-latest_filter-road/raw/slice-5_edges.geoparquet")
>>> df
                                                geometry  osm_way_id  segment_id  ... end_node_longitude end_node_latitude end_node_degree
0      LINESTRING (31.33793 30.09310, 31.33710 30.092...     4943484           0  ...          31.337104         30.092643               1
1      LINESTRING (31.33296 30.09032, 31.33294 30.090...     4943484           1  ...          31.332943         30.090311               1
2      LINESTRING (31.33055 30.08897, 31.33013 30.088...     4943484           2  ...          31.330126         30.088710               1
3      LINESTRING (31.39918 30.12437, 31.39922 30.124...     4943485           0  ...          31.399215         30.124412               1
4      LINESTRING (31.39960 30.12548, 31.40013 30.127...     4943485           1  ...          31.400126         30.127065               1
...                                                  ...         ...         ...  ...                ...               ...             ...
95549  LINESTRING (31.39697 30.07334, 31.39697 30.073...  1102620140           0  ...          31.396965         30.073378               1
95550  LINESTRING (31.39693 30.07344, 31.39692 30.073...  1102620140           1  ...          31.396918         30.073455               1
95551  LINESTRING (31.39664 30.07348, 31.39661 30.073...  1102620140           2  ...          31.396612         30.073454               1
95552  LINESTRING (31.39662 30.07321, 31.39663 30.073...  1102620140           3  ...          31.396631         30.073202               1
95553  LINESTRING (31.39685 30.07317, 31.39689 30.073...  1102620140           4  ...          31.396885         30.073191               1

[95554 rows x 16 columns]

>>> df[df.start_node_reference.isna()]
                                                geometry  osm_way_id  segment_id  ... end_node_longitude end_node_latitude end_node_degree
1359   LINESTRING (31.19696 30.35887, 31.19694 30.35905)    27649088           0  ...          31.196937         30.359048               2
1754   LINESTRING (31.23456 30.05621, 31.23484 30.056...    28851338           0  ...          31.234841         30.056289               1
1798   LINESTRING (31.26904 29.96933, 31.26928 29.969...    28917006           0  ...          31.269284         29.969280               1
2995   LINESTRING (31.29969 30.12808, 31.29968 30.128...    33622398           0  ...          31.299682         30.128017               1
3242   LINESTRING (32.86118 31.08451, 32.86212 31.083...    34423411           0  ...          32.862116         31.083364               1
...                                                  ...         ...         ...  ...                ...               ...             ...
95302  LINESTRING (31.30386 29.96250, 31.30381 29.962...  1094777357           0  ...          31.303808         29.962600               1
95491  LINESTRING (31.08164 29.34396, 31.08023 29.344...  1098484176           0  ...          31.080232         29.344363               1
95505  LINESTRING (31.23913 30.15031, 31.23902 30.150...  1099842244           0  ...          31.239017         30.150265               1
95516  LINESTRING (30.83173 31.23758, 30.83245 31.237...  1102285534           0  ...          30.832454         31.237803               1
95519  LINESTRING (29.85719 31.09891, 29.85702 31.10015)  1102358629           0  ...          29.857021         31.100152               2

[2028 rows x 16 columns]

>>> from shapely.geometry.point import Point
>>> import matplotlib.pyplot as plt
>>> mask = df.start_node_reference.isna()
>>> pts = [Point(x, y) for x, y in zip(df[mask].start_node_longitude, df[mask].start_node_latitude)]
>>> gdf = gpd.GeoDataFrame({"geometry": pts})
>>> ax = gdf.plot()
>>> ax.set_title("start_nodes with reference=-1")
>>> plt.show()
```

[Plot of node locations]('node_reference_problem.png')
