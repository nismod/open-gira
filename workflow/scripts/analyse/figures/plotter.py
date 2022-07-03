"""Plotter for caribbean region"""


import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os


try:
    output_dir = snakemake.params['output_dir']
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
    linewidth = .55
    output_dir = 'results'
    #raise RuntimeError("Please use snakemake to define inputs")

max_dist = .5


#print('Uncomment line below and world.plot (for testing speedup, was off)')
world_file = os.path.join(output_dir, 'input', 'adminboundaries', 'gadm36_levels.gpkg')
#world_file = """C:/Users/maxor/Documents/PYTHON/GIT/open-gira/results/input/adminboundaries/gadm36_levels.gpkg"""
world = gpd.read_file(world_file, layer=0)

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


# filter data
if linewidth:  # if not None
    data['len'] = [len(list(geom[1].coords)) for geom in data.geometry.items()]  # number of points on linestring
    data['geolen'] = [geom[1].length for geom in data.geometry.items()]
    data = data[~((data.len==2) & (data.geolen > max_dist))]  # filter overly long for visual



world.plot(ax=ax, color=(.9,.9,.9))  # TODO
data.plot(column=metric, ax=ax, linewidth=linewidth, legend=True, cmap=cmap, vmax=vmax, vmin=vmin, cax=cax, legend_kwds={'label':metric})
#world.boundary.plot(ax=ax, linewidth=.01, color=(0.5,0.5,0.5))
fig.set_figheight(11)
fig.set_figwidth(20)
plt.savefig(output, dpi=200, bbox_inches='tight')
print(f'Saved {output}')
