"""
Creates the network from the plants and targets data
"""

from importing_modules import *
from process_power_functions import *

#codes = COUNTRY_CODES
codes = sys.argv[1]
print(codes)
codes = ast.literal_eval(codes)  # convert to list


# codes = ["PHL"]
# changedir()
#%%
print("combining all intermediate files")
initial = False
for ii, code in enumerate(codes):
    print(code)
    plantsfile = os.path.join("data","processed",f"{code.upper()}_plants.csv")
    targetsfile = os.path.join("data","processed",f"{code.upper()}_targets.csv")

    plants_country = pd.read_csv(plantsfile)
    targets_country = pd.read_csv(targetsfile)

    if initial == False and len(plants_country) != 0:
        plants = plants_country
        targets = targets_country
        initial = True

    elif initial == True:
        print("appending")
        plants = plants.append(plants_country)
        targets = targets.append(targets_country)

plants['geometry'] = plants['geometry'].apply(wkt.loads)  # TODO wont need this if using gpkg (later fix, works for now)
plants = gpd.GeoDataFrame(plants)

targets['geometry'] = targets['geometry'].apply(wkt.loads)
targets = gpd.GeoDataFrame(targets)

#print("resetting indices")  # not needed
plants = plants.reset_index().copy()
plants['id'] = [f"source_{i}" for i in range(len(plants))]
plants = plants[['id', 'source_id', 'capacity_mw', 'type', 'geometry']]

targets = targets.reset_index().copy()
targets['id'] = [f"target_{i}" for i in range(len(targets))]


#%%
print("\n-- countries processed --\n")
# write to file

print("writing to world_plants.gpkg")
plants.to_file(os.path.join("data","processed","world_plants.gpkg"),driver='GPKG')
timer(start)

print("writing to world_targets.gpkg")
if len(targets) == 0:
    raise RuntimeError("targets file is empty, program will quit.")
targets.drop(columns=['centroid']).to_file(os.path.join("data","processed","world_targets.gpkg"), driver='GPKG')
timer(start)

print("processing targets GeoData")
targets = targets[['id', 'type', 'centroid']] \
    .rename(columns={'centroid': 'geometry'})
timer(start)

# Combine to nodes
print("combining sources and sinks")
nodes = plants.append(targets).reset_index(drop=True)  # sources and targets located
timer(start)

# Edges
print("getting gridfinder lines")
if len(codes) == 1:  # TODO for testing purposes, one country
    edges = get_lines(codes[0])
else:
    edges = get_lines()
timer(start)

print('processing lines')
edges['type'] = 'transmission'
edges['id'] = targets.reset_index()['index'].apply(lambda i: f"edge_{i}")

timer(start)
print("writing to world_edges.gpkg")
edges.to_file(
    os.path.join("data","processed","world_edges.gpkg"),
    driver='GPKG'
)
timer(start)

print("processing edges GeoData")
edges = edges[['id', 'source_id', 'type', 'geometry']]
timer(start)

# Process network
print("creating network")
network = snkit.network.Network(nodes, edges)
timer(start)

# fix str when should be Point(# #)
network.nodes['geometry'] = [sw.loads(x) if type(x) == str else x for x in network.nodes['geometry']]
network.edges['geometry'] = [sw.loads(x) if type(x) == str else x for x in network.edges['geometry']]



# Connect power plants
print("connecting powerplants")
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
    warnings.simplefilter(action='ignore', category=FutureWarning)

    print("processing network")
    network = snkit.network.split_multilinestrings(network)
    timer(start)

    geod = Geod(ellps="WGS84")
    edge_limit = 20_000 # meters

    network = snkit.network.link_nodes_to_nearest_edge(
        network,
        lambda node, edge: node.type == 'source' and geod.geometry_length(edge.geometry) < edge_limit)

network.nodes.loc[network.nodes.id.isnull(), 'type'] = 'conn_source'
timer(start)

print("resetting indices")
network.nodes['id'] = network.nodes.reset_index() \
    .apply(
        lambda row: f"conn_source_{row['index']}" if type(row.id) is float else row.id,
        axis=1
    )

# Connect targets
print("connecting targets")
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
    warnings.simplefilter(action='ignore', category=FutureWarning)

    network = snkit.network.link_nodes_to_nearest_edge(
        network,
        lambda node, edge: node.type == 'target' and geod.geometry_length(edge.geometry) < edge_limit)
network.nodes.loc[network.nodes.id.isnull(), 'type'] = 'conn_target'
timer(start)

print("resetting indices")
network.nodes['id'] = network.nodes.reset_index() \
    .apply(
        lambda row: f"conn_target_{row['index']}" if type(row.id) is float else row.id,
        axis=1
    )
timer(start)

# Add nodes at line endpoints
print("adding nodes at line endpoints")
network = snkit.network.add_endpoints(network)
timer(start)

print("processing network")
network.nodes.loc[network.nodes.id.isnull(), 'type'] = 'intermediate'
timer(start)

print("resetting indices")
network.nodes['id'] = network.nodes.reset_index() \
    .apply(
        lambda row: f"intermediate_{row['index']}" if type(row.id) is float else row.id,
        axis=1
    )
timer(start)

# add from/to ids
print("adding from/to ids")
with warnings.catch_warnings():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    network = snkit.network.add_topology(network, id_col='id')
timer(start)

# output
out_fname = os.path.join("data","processed","world_network.gpkg")
print("writing edges to ", out_fname)
network.edges.to_file(out_fname, layer='edges', driver='GPKG')
timer(start)

print("writing nodes to ", out_fname)
network.nodes.to_file(out_fname, layer='nodes', driver='GPKG')
timer(start)

print("\nnote that ShapelyDepreciationWarnings (shapely) and FutureWarning (geopandas) were supressed. Also newer pandas versions can cause issues (1.1.5 works).")  # remove if fixed
end = time.time()
print(f"\nTime to run file: {(end - start)/60:0.2f} minutes")


