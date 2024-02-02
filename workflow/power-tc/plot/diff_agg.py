"""Takes in data/ file startswith1 minus startwith2

Best to have 1 as future

targets: eff pop ae and code
recon agg: recon cost ae and code
recon agg norm: recon cost ae frac norm and code
recon indiv: recon cost ae and link
"""
import os
import time

import geopandas as gpd

try:
    metric = snakemake.params["metric"]  # type: ignore
    merge_key = snakemake.params["merge_key"]  # type: ignore
    inputs = snakemake.input  # type: ignore
    output = snakemake.output  # type: ignore
except:
    raise RuntimeError("Please use snakemake to define inputs")

# first input is future
assert len(inputs) == 2
assert type(output) != list
output = str(output)

d = dict()  # master dict
g = dict()  # geom
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

print(f"len(d) pre = {len(d)}")
d = {k: v for k, v in d.items() if len(v) == 2}
print(f"len(d) post = {len(d)}")

p = dict()  # percentage
for k, v in d.items():
    d_orig = d[k][1]
    d[k] = d[k][0] - d[k][1]
    if d[k] != 0:
        p[k] = 100 * d[k] / d_orig
    else:
        p[k] = 0


gdf = gpd.GeoDataFrame({merge_key: g.keys(), "geometry": g.values()})
rem = "REMOVE"
gdf[metric] = gdf[merge_key].map(d).fillna(rem)
perc_metric = f"perc_{metric}"
gdf[perc_metric] = gdf[merge_key].map(p).fillna(rem)
b = [False if gdf.iloc[jj][metric] == rem else True for jj in range(len(gdf))]
gdf = gdf[b]

for col in gdf.columns:
    if col not in ["geometry", merge_key]:
        gdf[col] = gdf[col].astype(float)

output_folder = os.path.dirname(output)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

gdf.to_file(output, driver="GPKG")

print("done")
time.sleep(1)
