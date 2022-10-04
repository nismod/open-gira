"""Takes a gpkg file and aggregates to chosen level

Outputs a gpkg file with metrics as target features (can use QGIS to plot)
"""
import os

import fiona
import geopandas as gpd
from shapely.geometry import shape
from tqdm import tqdm

from .common_functions import ae, avg, sm

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
    metrics_target = snakemake.params["metrics_target"]  # type: ignore
    top_select = snakemake.params[  # type: ignore
        "top_select"
    ]  # select (in percent). Set to 100 for all
    increased_severity_sort = snakemake.params["increased_severity_sort"]  # type: ignore
    layer_num = snakemake.params["aggregate_level"]  # type: ignore
except:  # for testing only
    output_dir = "results"
    metrics_target = [
        "population_without_power",
        "effective_population",
        "affected_population",
        "mw_loss_storm",
        "f_value",
        "gdp_damage",
    ]
    top_select = 100
    increased_severity_sort = True
    layer_num = 1
    raise RuntimeError("Please use snakemake to define inputs")

increased_severity_sort_bool = str(increased_severity_sort)[0]  # either T or F


metrics_target_nof = metrics_target.copy()
metrics_target_nof.remove("f_value")
metrics_target_avg = [avg(metric) for metric in metrics_target]
metrics_target_sum = [
    sm(metric) for metric in metrics_target_nof
]  # f_value does not make sense for sum
metrics_target_ae = [
    ae(metric) for metric in metrics_target_nof
]  # f_value does not make sense for ae
metric_keys = metrics_target_sum + metrics_target_avg + metrics_target_ae


folder_agg = os.path.join(output_dir, "power_output", "statistics", "aggregate")
examine_file = os.path.join(
    folder_agg, f"targets_geo_top{top_select}{increased_severity_sort_bool}percent.gpkg"
)
quantile_file = gpd.read_file(examine_file, keep_default_na=False)
quantile_file.crs = "EPSG:4326"
countries_relevant = set(
    quantile_file.country.unique()
)  # set of countries in which to check (others irrelevant, significant computational improvement)

print(f"loading level {layer_num} data")
with fiona.open(
    os.path.join(output_dir, "input", "admin-boundaries", f"gadm36_levels.gpkg"),
    "r",
    layer=layer_num,
) as src_code:
    code_geoms = []
    code_GIDs = []
    for feature in src_code:
        if (
            feature["properties"]["GID_0"] in countries_relevant
        ):  # only include search in countries that contain targets
            code_geoms.append(shape(feature["geometry"]))
            code_GIDs.append(feature["properties"][f"GID_{layer_num}"])
    print("creating dataframe")
    code_geoms_gpd = gpd.GeoDataFrame(
        {"geometry": code_geoms, "code": code_GIDs}, crs="EPSG:4326"
    )

quantile_file["centre"] = quantile_file["geometry"].centroid


# code_geoms_gpd in quantile box
qf_bounds = quantile_file.bounds
qf_maxx = qf_bounds["maxx"].max()
qf_minx = qf_bounds["minx"].min()
qf_maxy = qf_bounds["maxy"].max()
qf_miny = qf_bounds["miny"].min()

cg_bounds = code_geoms_gpd.bounds
cg_maxx_bool = (
    cg_bounds["maxx"] < qf_minx
)  # True if max region x is less than min quantile x
cg_minx_bool = (
    cg_bounds["minx"] > qf_maxx
)  # True if min region x is more than max quantile x
cg_maxy_bool = (
    cg_bounds["maxy"] < qf_miny
)  # True if max region is less than min quantile y
cg_miny_bool = (
    cg_bounds["miny"] > qf_maxy
)  # True if min region is more than max quantile y
bool_keep = [
    False
    if cg_maxx_bool[x] + cg_minx_bool[x] + cg_maxy_bool[x] + cg_miny_bool[x] >= 1
    else True
    for x in range(len(cg_bounds))
]  # if 1 or higher, then there is at least one True element
code_geoms_gpd = code_geoms_gpd[bool_keep]

map_dict = {}

for geom_area in tqdm(
    code_geoms_gpd.itertuples(), total=len(code_geoms_gpd), desc="geom_intersect"
):

    bool_list = [
        True if geom_area.geometry.covers(c) else False for c in quantile_file.centre
    ]  # check contained (tehcnically covered as can lie on boundary, later is removed, so no issues with double counting) in geometry layer

    overlap_quantile = quantile_file[bool_list]
    if len(overlap_quantile) >= 1:
        remove_id = overlap_quantile["index"]
        if geom_area.code not in map_dict.keys():
            map_dict[geom_area.code] = {
                metric: [] for metric in metric_keys
            }  # dict(zip(metric_keys, [[]]*len(metric_keys)))
        for metric in metrics_target:
            if "f_value" in metric:
                map_dict[geom_area.code][avg(metric)] = overlap_quantile[
                    avg(metric)
                ].mean()
            else:
                map_dict[geom_area.code][sm(metric)] = overlap_quantile[
                    sm(metric)
                ].sum()
                map_dict[geom_area.code][avg(metric)] = overlap_quantile[
                    avg(metric)
                ].mean()
                map_dict[geom_area.code][ae(metric)] = overlap_quantile[
                    ae(metric)
                ].sum()

                if "population" in metric:
                    new_key = f"{metric}_annually-expected_region_fraction"
                    if new_key not in metric_keys:
                        metric_keys.append(new_key)
                    map_dict[geom_area.code][new_key] = (
                        overlap_quantile[ae(metric)].sum()
                        / overlap_quantile.population.sum()
                    )  # fraction of total target populations in the quantile file overlap (ie all targets in the geometry area)

        quantile_file = quantile_file[
            ~quantile_file["index"].isin(remove_id)
        ]  # remove classified targets


for metric in metric_keys:
    map_dict_indiv = {
        k: v[metric] for k, v in map_dict.items()
    }  # {code1: metric_value1, code2: metric_value2 ... }
    code_geoms_gpd[metric] = code_geoms_gpd["code"].map(map_dict_indiv).fillna(0)


new_file_name = examine_file.replace("percent", f"percent_aggregated_region")
if len(code_geoms_gpd) == 0:
    code_geoms_gpd = gpd.GeoDataFrame({"geometry": [None]})  # add dummy for snakemake
code_geoms_gpd.to_file(new_file_name, driver="GPKG")
print("top_select written to file")
