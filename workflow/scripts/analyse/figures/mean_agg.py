"""Takes anything in folder data/ and aggregated to admin keys"""


import os
import geopandas as gpd

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
    metric = snakemake.params["metric"]  # type: ignore
    merge_key = snakemake.params["merge_key"]  # type: ignore
    inputs = snakemake.input  # type: ignore
    output = snakemake.output  # type: ignore
except:
    raise RuntimeError("Please use snakemake to define inputs")


plot_path = os.path.join(output_dir, "power_figures")
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


assert type(output) != list
output = str(output)
d = dict()
g = dict()
for ii, file in enumerate(inputs):

    mean_file = gpd.read_file(file)
    for jj in range(len(mean_file)):
        code = mean_file.iloc[jj][merge_key]
        metric_val = mean_file.iloc[jj][metric]
        if code not in d:
            d[code] = [metric_val]
            g[code] = mean_file.iloc[jj]["geometry"]
        else:
            d[code] = d[code] + [metric_val]


for k, v in d.items():
    d[k] = sum(d[k]) / len(d[k])

gdf = gpd.GeoDataFrame({merge_key: g.keys(), "geometry": g.values()})
gdf[metric] = gdf[merge_key].map(d).fillna(0)

output_folder = os.path.dirname(output)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

gdf.to_file(output, driver="GPKG")

print("done")
time.sleep(1)
