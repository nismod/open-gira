"""Plotter for caribbean region"""


import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    plotfile = str(snakemake.input)
    output = str(snakemake.output)
    metric = snakemake.params['metric']
    vmax = snakemake.params['vmax']
    vmin = snakemake.params['vmin']
    cmap = snakemake.params['cmap']
    linewidth = snakemake.params['linewidth']
except:
    plotfile = "C:/Users/maxor/Documents/PYTHON/GIT/open-gira/results/power_figures/intermediate_files/difference_link_reconstruction_cost_annual_expected.gpkg"
    output = "C:/Users/maxor/Documents/PYTHON/GIT/open-gira/results/power_figures/fig1.png"
    metric = 'effective_population_anually-expected'
    metric = 'perc_reconstruction_cost_annual_expected'
    vmax = 30
    vmin = -30
    cmap = 'RdBu_r'
    linewidth = 1
    raise RuntimeError("Please use snakemake to define inputs")


print('Uncomment line below and world.plot (for testing speedup, was off)')
world = gpd.read_file("C:/Users/maxor/Documents/PYTHON/GIT/open-gira/results/input/adminboundaries/gadm36_levels.gpkg", layer=0)

minx = -87
maxx = -63
maxy = 24
miny = 16




fig, ax = plt.subplots(1, 1, frameon = False)

data = gpd.read_file(plotfile)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)


ax.set_xlim(minx, maxx)
ax.set_ylim(miny, maxy)

world.plot(ax=ax, color=(.9,.9,.9))  # TODO
data.plot(column=metric, ax=ax, linewidth=linewidth, legend=True, cmap=cmap, vmax=vmax, vmin=vmin, cax=cax, legend_kwds={'label':metric})
fig.set_figheight(11)
fig.set_figwidth(20)
plt.savefig(output, dpi=200)
print(f'Saved {output}')
