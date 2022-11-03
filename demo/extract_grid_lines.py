import geopandas

df = geopandas.read_file("results/input/gridfinder/grid.gpkg")
df.reset_index(inplace=True, names="edge_id")
df.edge_id = df.edge_id.apply(lambda e: f"grid_{e}")
df.to_parquet("results/grid.geoparquet")

# python workflow/scripts/extract_grid_lines.py
# python workflow/scripts/intersection.py results/grid.geoparquet "results/input/storm-ibtracs/fixed/*" results/grid_split.geoparquet

# set up damages curve and cost estimate
# run return period damages
# run expected damages
# aggregate back to line
# aggregate to regions
