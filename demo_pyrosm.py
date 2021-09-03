#!/usr/bin/env python
# coding: utf-8

# In[37]:


import os
from glob import glob


import geopandas
import pandas
import pyrosm
import rasterio
import snail

from snail.intersections import split
from snail.intersections import get_cell_indices
from tqdm.notebook import tqdm


# In[2]:


adm1_name = 'bangladesh'


# In[8]:


data_folder = '/tmp/mert2014'


# In[3]:


osm = pyrosm.OSM(os.path.join(data_folder,'osm',f'{adm1_name}-latest-highway.osm.pbf'))


# In[5]:


nodes, edges = osm.get_network(nodes=True, network_type="driving")


# In[12]:


core = (
    'motorway_link',
    'motorway',
    'trunk_link',
    'trunk',
    'primary_link',
    'primary',
    'secondary_link',
    'secondary',
    'tertiary_link',
    'tertiary',
)
core_edges = edges[edges.highway.isin(core)]


# In[13]:


len(core_edges), len(edges)


# In[14]:


select_columns = [
    'bridge', 'highway', 'lanes', 'maxspeed', 'oneway',
    'smoothness', 'surface', 'tracktype', 'tunnel', 'width', 
    'id', 'name', 'osm_type', 'geometry', 'u', 'v', 'length'
]
core_edges = core_edges[select_columns]


# In[18]:


get_ipython().run_cell_magic('timeit', '', "core_edges.to_file(os.path.join(data_folder, 'osm', f'{adm1_name}-roads-core.gpkg'), driver='GPKG')")


# In[15]:


get_ipython().run_cell_magic('timeit', '', "core_edges.to_file(os.path.join(data_folder, 'osm', f'{adm1_name}-roads-core.fgb'), driver='FlatGeobuf')")


# In[16]:


get_ipython().run_cell_magic('timeit', '', "core_edges.to_parquet(os.path.join(data_folder, 'osm', f'{adm1_name}-roads-core.geoparquet'))")


# In[17]:


get_ipython().run_cell_magic('timeit', '', "core_edges.to_feather(os.path.join(data_folder, 'osm', f'{adm1_name}-roads-core.geofeather'))")


# In[19]:


# Write direct from pyrosm driving
#edges.to_file(os.path.join(data_folder, 'osm', f'{adm1_name}-roads.gpkg'), driver='GPKG', layer='edges')
#nodes.to_file(os.path.join(data_folder, 'osm', f'{adm1_name}-roads.gpkg'), driver='GPKG', layer='nodes')


# In[20]:


# Read from file written above
#core_edges = geopandas.read_file(os.path.join('data', 'osm', f'{adm1_name}-roads-core.gpkg'))


# In[21]:


raster_data = rasterio.open(os.path.join('..', 'aqueduct', 'inuncoast_historical_nosub_hist_rp0050_0.tif'))


# In[97]:


get_ipython().run_cell_magic('timeit', '', "core_splits = []\nfor edge in tqdm(core_edges.itertuples()):\n    splits = split(\n        edge.geometry,\n        raster_data.width,\n        raster_data.height,\n        list(raster_data.transform),\n    )\n    for s in splits:\n        core_splits.append({\n            'id': edge.id,\n            'geometry': s\n        })\ncore_splits = geopandas.GeoDataFrame(core_splits)")


# In[98]:


len(core_edges), len(core_splits)


# In[23]:


core_edges[['id','geometry']].head(50).tail()


# In[24]:


tqdm.pandas()


# In[96]:


get_ipython().run_cell_magic('timeit', '', "core_splits['cell_index'] = core_splits.geometry.apply(\n    lambda geom: list(get_cell_indices(geom, raster_data.width, raster_data.height, list(raster_data.transform))))")


# In[27]:


core_splits.head()


# In[28]:


band = raster_data.read(1)


# In[29]:


get_ipython().run_cell_magic('timeit', '', "core_splits['inuncoast_historical_nosub_hist_rp0050_0'] = core_splits.cell_index.apply(lambda i: band[i[1], i[0]])")


# In[30]:


core_splits


# In[31]:


fnames = glob('../aqueduct/*.tif')
fnames[0]


# In[52]:


coastal = []
river = []
for fname in fnames:
    fname = os.path.basename(fname)
    colname = fname[:-4]
    if 'coast' in colname:
        # inuncoast_{climatescenario}_{subsidence}_{year}_{returnperiod}_{projection}.tif
        try:
            _, clim, sub, y, rp, proj = colname.split("_")
        except ValueError:
            _, clim, sub, y, rp, _, _, proj = colname.split("_")
        if proj == "0":
            proj = "95"
        if y == "hist":
            y = 2010
        coastal.append({
            "key": colname,
            "climate_scenario": clim,
            "subsidence": sub,
            "year": int(y),
            "return_period": int(rp[2:]),
            "sea_level_rise_percentile": int(proj),
            "filename": fname,
        })
    else:
        # inunriver_{climatescenario}_{model}_{year}_{returnperiod}.tif
        _, clim, model, y, rp = colname.split("_")
        if y == "hist":
            y = 2010
        
        river.append({
            "key": colname,
            "climate_scenario": clim,
            "model": model.replace("0",""),
            "year": int(y),
            "return_period": int(rp[2:]),
            "filename": fname,
        })
coastal = pandas.DataFrame(coastal)
river = pandas.DataFrame(river)


# In[53]:


coastal.to_csv('aqueduct_coastal.csv')


# In[54]:


river.to_csv('aqueduct_river.csv')


# In[60]:


def associate_raster(df, key, fname, band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df.cell_index.apply(lambda i: band_data[i[1], i[0]])


# In[57]:


river.year.value_counts()


# In[59]:


subset = river[river.year.isin((1980, 2080)) & river.return_period.isin((50, 100, 500, 1000))]
len(subset)


# In[91]:


get_ipython().run_cell_magic('timeit', '', "associate_raster(core_splits, 'inunriver_rcp8p5_00IPSL-CM5A-LR_2080_rp00050', os.path.join(data_folder, 'aqueduct', 'inunriver_rcp8p5_00IPSL-CM5A-LR_2080_rp00050.tif'))")


# In[92]:


# do I/O stuff outside of timeit loop
dataset = rasterio.open(os.path.join(data_folder, 'aqueduct', 'inunriver_rcp8p5_00IPSL-CM5A-LR_2080_rp00050.tif'))
band_data = dataset.read(1)


# In[95]:


len(core_splits)


# In[94]:


get_ipython().run_cell_magic('timeit', '', "core_splits['inunriver_rcp8p5_00IPSL-CM5A-LR_2080_rp00050'] = core_splits.cell_index.apply(lambda i: band_data[i[1], i[0]])")


# In[86]:


for raster in subset.itertuples():
    associate_raster(core_splits, raster.key, os.path.join(data_folder, 'aqueduct', raster.filename))


# In[87]:


get_ipython().run_cell_magic('timeit', '', "core_splits.drop(columns='geometry').to_csv(os.path.join(data_folder, 'outputs', 'core_splits.csv.gz'))")


# In[81]:


get_ipython().run_cell_magic('timeit', '', "pandas.DataFrame(core_splits.drop(columns=['geometry'])) \\\n    .to_parquet(os.path.join(data_folder, 'outputs', 'core_splits.parquet'))")


# In[88]:


get_ipython().run_cell_magic('timeit', '', "core_splits[['id','geometry']].to_file(os.path.join(data_folder,  'outputs', 'core_splits.gpkg'), driver='GPKG')")


# In[89]:


get_ipython().run_cell_magic('timeit', '', "core_splits.to_parquet(os.path.join(data_folder,  'outputs', 'core_splits.geoparquet'))")


# In[90]:


get_ipython().run_cell_magic('timeit', '', "core_splits.to_feather(os.path.join(data_folder,  'outputs', 'core_splits.geofeather'))")

