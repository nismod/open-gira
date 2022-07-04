"""Plotter for caribbean region"""


import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from shapely.geometry import LineString

try:
    output_dir = snakemake.params["output_dir"]
    plotfile = str(snakemake.input)
    output = str(snakemake.output)
    metric = snakemake.params["metric"]
    vmax = snakemake.params["vmax"]
    vmin = snakemake.params["vmin"]
    cmap = snakemake.params["cmap"]
    linewidth = snakemake.params["linewidth"]
    legend_name = snakemake.params["legend_name"]
except:
    plotfile = "N/A"
    output = "N/A"
    metric = "effective_population_anually-expected"
    metric = "perc_reconstruction_cost_annual_expected"
    vmax = 30
    vmin = -30
    cmap = "RdBu_r"
    linewidth = 0.55
    output_dir = "results"
    legend_name = "LEG NAME"
    raise RuntimeError("Please use snakemake to define inputs")

max_dist = 0.5


world_file = os.path.join(output_dir, "input", "adminboundaries", "gadm36_levels.gpkg")
world = gpd.read_file(world_file, layer=0)

minx = -87
maxx = -63
maxy = 24
miny = 16


def eval_dist_lst(linestring_df):
    """"Evaluate the coordinate dataframe and returns distance list in km"""

    lst = []
    for ii in range(len(linestring_df)):
        dist_tot = 0
        if type(linestring_df.iloc[ii].geometry) == type(
            LineString([[1, 2], [3, 4]])
        ):  # check is linestring
            line_coords = list(
                linestring_df.iloc[ii].geometry.coords
            )  # extract the coordinates of a row
            dist_tot = len(line_coords)

        else:  # multistring
            for ms in range(len(linestring_df.iloc[ii]["geometry"].geoms)):
                line_coords = list(
                    linestring_df.iloc[ii].geometry[ms].coords
                )  # extract the coordinates of a row
                dist_tot += len(line_coords)

        lst.append(dist_tot)

    return lst


fig, ax = plt.subplots(1, 1, frameon=False)

data = gpd.read_file(plotfile)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)


ax.set_xlim(minx, maxx)
ax.set_ylim(miny, maxy)


# filter data
if linewidth:  # if not None
    data["len"] = eval_dist_lst(data)  # number of points on linestring
    data["geolen"] = [geom[1].length for geom in data.geometry.items()]
    data = data[
        ~((data.len == 2) & (data.geolen > max_dist))
    ]  # filter overly long for visual


if (
    vmax == "evenmax" and vmin == "evenmin"
):  # ensure either side of zero is equal distance from vmax and vmin
    vmax = max(data[metric].max(), -data[metric].min())
    vmin = max(-data[metric].max(), data[metric].min())


if vmax == "max":
    vmax = data[metric].max()
if vmin == "min":
    vmin = data[metric].min()


world.plot(ax=ax, color=(0.9, 0.9, 0.9))
data.plot(
    column=metric,
    ax=ax,
    linewidth=linewidth,
    legend=True,
    cmap=cmap,
    vmax=vmax,
    vmin=vmin,
    cax=cax,
    legend_kwds={"label": legend_name},
)
data.boundary.plot(ax=ax, linewidth=0.1, color=(0.5, 0.5, 0.5))
fig.set_figheight(11)
fig.set_figwidth(20)
plt.savefig(output, dpi=200, bbox_inches="tight")
print(f"Saved {output}")
