"""Calculates the gdp losses for a given storm
"""

import os
import pandas as pd
import geopandas as gpd
import time
from datetime import date
from shapely.geometry import LineString
import json
import sys
from damage_calculator import applythreshold
import networkx as nx
import numpy as np
import shapely.wkt as sw





if 'linux' not in sys.platform:  # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

try:
    region = snakemake.params["region"]
    sample = snakemake.params["sample"]
    nh = snakemake.params["nh"]
except:
    print('RUNNING ON FIXED INPUTS')
    region = 'NA'
    sample = '0'
    nh = '0_11_11'

def isNone(df):
    """Checks if dataframe contains solely the None row (required for snakemake and gpkg files)"""
    if len(df) == 1 and df['id'].iloc[0] == None:
        return True
    else:
        return False


def read_edges_make_unique(fname):
    """Read edges and add unique link column"""
    edges = gpd.read_file(fname, layer="edges")

    if len(edges) == 1:
        if edges["id"].iloc[0] == None:  # no edges
            edges["link"] = None
    else:
        # create unique link id from from/to node ids
        edges["link"] = edges.apply(
            lambda e: "__".join(sorted([e.from_id, e.to_id])), axis=1
        )
    return edges


def read_network(fname):
    """Read network and add relevant attributes"""
    # read edges
    edges = read_edges_make_unique(fname)

    # read nodes
    nodes = gpd.read_file(fname, layer="nodes")

    # add gdp to target nodes
    targets = gpd.read_file(fname.replace("network", "targets"))
    nodes = nodes.merge(targets[["id", "gdp"]], on="id", how="left")

    # add capacity to plant/generation/source nodes
    plants = gpd.read_file(fname.replace("network", "plants"))
    if nodes["capacity_mw"].all() == plants["capacity_mw"].all():
        nodes = nodes.drop(columns="capacity_mw")  # prevent merge issue
    nodes = nodes.merge(plants[["id", "capacity_mw"]], on="id", how="left")

    return nodes, edges


def create_graph(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from(
        (n.id, {"id": n.id, "type": n.type, "gdp": n.gdp, "capacity_mw": n.capacity_mw})
        for n in nodes.itertuples()
    )
    G.add_edges_from(
        (e.from_id, e.to_id, {"id": e.link})
        for e in edges.itertuples()
    )
    return G


def box_connectors(box_id):
    """Loads and returns the connector dictionary"""
    with open(
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            f"{box_id}",
            f"connector_{box_id}.txt",
        ),
        "r",
    ) as file_ex:
        connector_adj = json.load(file_ex)
    return connector_adj

def network_name(box_id):
    fname = os.path.join(
        "data", "processed", "all_boxes", box_id, f"network_{box_id}.gpkg"
    )
    return fname

def component_select(box_id, node_set, network_dict):
    """Finds the connected component(s) in box_id which contain(s) node_set
    Input:
        box_id: str
        node_set: set of nodes which are required to be in the connected component(s)
        network_dict: dictionary of (nodes, edges) corresponding to the box_id keys (to reduce time on loading, only do once)
    Output:
        nodes: dataframe of nodes connected to node_set
        edges: dataframe of edges connected to edge_set
        node_set: updated set of nodes (does not anymore include the nodes which have been found in nodes, edges (above))
        conn_boxes: set of boxes which need to be examined still

    """
    connectors = box_connectors(box_id)
    conn_boxes = set()
    nodes = gpd.GeoDataFrame()
    edges = gpd.GeoDataFrame()


    #print(box_id)

    ## Option 1
    #nodes_init, edges_init = read_network(network_name(box_id))  # load the network for the boxdamaged component

    ## Option 2
    nodes_init, edges_init = network_dict[box_id]


    if isNone(nodes_init) or isNone(edges_init):  # if None row, return empty
        return nodes, edges, node_set, conn_boxes


    G_init = create_graph(nodes_init, edges_init)  # create a graph
    components_init = list(nx.connected_components(G_init))  # split into components
    for comp_init in components_init:  # for each component
        if any(x in comp_init for x in node_set):  # if any node from the node_set and current node list is in this component

            nodes_comp = nodes_init[nodes_init['id'].isin(comp_init)]  # nodes in this component
            edges_comp = edges_init[
                np.maximum(
                    edges_init["from_id"].isin(comp_init),
                    edges_init["to_id"].isin(comp_init),
                )]  # edges in this component

            nodes = nodes.append(nodes_comp)  # keep the nodes of this component
            edges = edges.append(edges_comp)  # keep only edges of this component (note max(True, False) = True)

            ## note which nodes have been examined and which remain ##
            node_set_in_comp_init = {node for node in node_set if node in comp_init}  # nodes in component examined now, so will be removed
            conn_dict = {k: [v_conn['to_id'] for v_conn in v_lst if v_conn['from_id'] in comp_init.difference(node_set_in_comp_init)] for k, v_lst in connectors.items()}  # add the nodes (to_id) if the from_id is in the component (but not if connected to the previous box)
            node_set = node_set.difference(node_set_in_comp_init)  # node_set now contains all outer search nodes still to be examined
            conn_boxes = set(conn_dict.keys())  # note the boxes to still be examined
            conn_new = set().union(*conn_dict.values())  # note the nodes still to be examined
            node_set = node_set.union(conn_new)  # update to include new node connections


            ## add the links for the examined cross-box connections ##
            cross_box_links = []  # list of dictionary with the relevant cross-box connections
            for connectors_indiv in connectors.values():
                cross_box_links += [conn_link_dict for conn_link_dict in connectors_indiv if conn_link_dict['from_id'] in node_set_in_comp_init]  # include the links where the 'from_id' is in the nodes in comp_init. That means the link travels into the box currently being examined
            if len(cross_box_links) != 0:
                update_gdf = gpd.GeoDataFrame(cross_box_links)
                update_gdf['geometry'] = update_gdf['geometry'].apply(sw.loads)  # make geometry type
                edges = edges.append(update_gdf)  # add the new links

    return nodes, edges, node_set, conn_boxes







def combine_networks(
    edge_damaged,
):
    """finds adjacent boxes which are connected to original box_id and adds to nodes and edges.
    Returns all nodes and edges (both connected into base_box). edge_damaged is dataframe input containing link and box_id and id"""


    box_id_orig = edge_damaged.box_id
    fname = os.path.join(
        "data", "processed", "all_boxes", box_id_orig, f"network_{box_id_orig}.gpkg"
    )
    nodes_orig, edges_orig = read_network(fname)  # load the network for the box which contains the damaged component
    G_orig = create_graph(nodes_orig, edges_orig)  # create a graph
    components_orig = list(nx.connected_components(G_orig))  # split into components
    for comp_orig in components_orig:  # for each component
        if any(x in comp_orig for x in [edge_damaged.from_id, edge_damaged.to_id]):  # if either end (from_id or to_id) of the damaged component is in this component


            nodes = nodes_orig[nodes_orig['id'].isin(comp_orig)]  # keep only the nodes of this component
            edges = edges_orig[
                np.maximum(
                    edges_orig["from_id"].isin(comp_orig),
                    edges_orig["to_id"].isin(comp_orig),
                )]  # keep only edges of this component (note max(True, False) = True)



            break  # edge_damaged.id will not be in remaining components, stop looking



    ## For the original box, note the connectors to conn_set ##
    conn_dict = {key_box: {conn_dict['to_id'] for conn_dict in conn_lst if conn_dict['from_id'] in nodes.id.values} for key_box, conn_lst in box_connectors(box_id_orig).items()}  # if any link starts in the current nodes df, then add it to the dictionary
    conn_set = set().union(*conn_dict.values())  # conn set will contain a set of links which have to be explored

    to_examine = set(conn_dict.keys())  # set of boxes to examine

    c = 0


    ## comment out below for if better memory available)
    network_dict = dict()  # dictionary {box_id1: (nodes_of_box_id1, edges_od_box_id1), box_id2, ...}
    while len(conn_set) >= 1:  # while there are still links to be examined



        box_id_examine = min(to_examine)  # pick one from to_examine (which is irrelevant)


        if box_id_examine not in network_dict.keys():  # if not in the dictionary, then add it
            network_dict[box_id_examine] = read_network(network_name(box_id_examine))
            print('added to dict')
        else:
            print('already in dict')



        c += 1
        if c%25 == 0:
            print(f"Examined {c} boxes")
        to_examine = to_examine.difference({box_id_examine})

        #print('examining ', box_id_examine)
        nodes_examine, edges_examine, conn_set, to_examine_newboxes = component_select(box_id_examine, conn_set, network_dict)  # extract the network components which connect to any node in conn_set, thenupdate conn_set to not include these are more
        nodes = nodes.append(nodes_examine)  # add
        edges = edges.append(edges_examine)  # add
        to_examine = to_examine.union(to_examine_newboxes)  # update
        #print(len(conn_set))


        # # TODO for testing
        # c_file = os.path.join('tester', 'expansion', f'expansion_{edge_damaged.id}_{c}.gpkg')
        # edges.to_file(c_file, layer="edges", driver="GPKG")
        # nodes.to_file(c_file, layer="nodes", driver="GPKG")


    return nodes, edges


def target_mapper(feature, targets, nodes):
    """Maps feature (of nodes df) to targets df. Will remove rows in which map has no value"""
    target_mapper_dict = {target_id: feature for target_id, feature, node_type in zip(nodes["id"].values, nodes[feature].values, nodes["type"].values) if node_type=='target'}  # dictionary {target1: f_value_of_target1, target2: f_value_of_target2, ... }
    targets[feature] = targets['id'].map(target_mapper_dict).fillna('remove_me')  # map to targets, note which are not connected (with 'remove_me')
    return targets[targets.f_value!='remove_me']  # remove unwanted targets



# set iteration variables
routeid_damaged = set() # marks which source_sink routes already damaged -> do not double count. Form: {(source1, target1), (source2, target2), ...}
totdamage = 0  # total damage (ensuring no double counts)
targetsdamaged = {}  # if the operation value of the target is desired
edges_affected = gpd.GeoDataFrame()
targets = gpd.GeoDataFrame()
polys_affected = gpd.GeoDataFrame()



print(f"{nh}: loading data")
# print('loading tracks')
stormfile = os.path.join(
    "data",
    "stormtracks",
    "events",
    f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
)
TC = pd.read_csv(stormfile, header=None)
TC.columns = [
    "year",
    "month",
    "number",
    "step",
    "basin",
    "lat",
    "lon",
    "pressure",
    "wind",
    "radius",
    "cat",
    "landfall",
    "dis_land",
]  # https://www.nature.com/articles/s41597-020-0381-2.pdf
TC = TC[["year", "number", "lat", "lon"]]
TC["number_hur"] = (
    str(sample)
    + "_"
    + TC["year"].astype(int).astype(str)
    + "_"
    + TC["number"].astype(int).astype(str)
)

# print("Loading wind")
windfile = os.path.join(
    "data",
    "intersection",
    "storm_data",
    "all_winds",
    region,
    sample,
    f"TC_r{region}_s{sample}_n{nh}.csv",
)
if not os.path.isfile(windfile):
    raise OSError(
        f"Wind file should exist but doesnt (TC_r{region}_s{sample}_n{nh}.csv)"
    )

winds_ev_all = pd.read_csv(windfile)
assert len(winds_ev_all) != 0
if not isNone(windfile):




    print(f"Investigating storm {nh}")
    winds_ev_filtered = applythreshold(winds_ev_all)
    ID_affected = list(winds_ev_filtered["ID_point"])
    box_id_affected = winds_ev_filtered["box_id"].unique()

    # print("- grid")
    grid_data = gpd.read_file(
        os.path.join("data", "intersection", "regions", f"{region}_unit.gpkg")
    )

    polys_affected = grid_data[grid_data["ID_point"].isin(ID_affected)].rename(columns={'box_id':'box_id_poly'})


    start = time.time()


    for jj, box_id in enumerate(box_id_affected):  # extract the damaged edges using this for loop
        print(f"-- Inspecting for damage {jj+1}/{len(box_id_affected)} -- {box_id}")

        box_edges = gpd.read_file(
            os.path.join(
                "data", "processed", "all_boxes", box_id, f"network_{box_id}.gpkg"
            ),
            layer="edges",
        )

        box_edges_affected = box_edges.overlay(
            polys_affected, how="intersection"
        )  # keeps edges that are affected grid points (only a part has to be in)


        box_edges_affected["link"] = box_edges_affected.apply(lambda e: "__".join(sorted([e.from_id, e.to_id])), axis=1)  # consistent naming

        edges_affected = edges_affected.append(box_edges_affected)  # add to master list of damaged edges


    edges = gpd.GeoDataFrame(columns=['link'])
    nodes = gpd.GeoDataFrame()

    print("Starting network connection expansion")
    startx = time.time()
    for edge_damaged in edges_affected.itertuples():


        if edge_damaged.link not in set(edges['link']):
            print('New link')
            #print('Searching Network')
            nodes_new, edges_new = combine_networks(edge_damaged)
            s1 = time.time()
            nodes = nodes.append(nodes_new)
            s2 = time.time()
            #print(f"Time for s1 is {s2 - s1}")
            edges = edges.append(edges_new)
            #print(f"Time for s2 is {time.time() - s2}")

    print(f'Search took {round((time.time() - startx)/60,1)} mins')


    ## Now all the damaged edges can be found in nodes & edges
    G = create_graph(nodes, edges)

    components = list(nx.connected_components(G))



    ## Set nominal values ##
    nodes['nominal_mw'] = 0  # set base to zero
    nominal_dict = dict()  # dictionary: {component1: {'nominal_mw': nominal_mw_1, 'nominal_gdp': nominal_gdp_1}, ... }
    for ii, component in enumerate(components):
        component_nodes = nodes[nodes['id'].isin(component)]
        nodes.loc[nodes.id.isin(component), 'component'] = ii
        total_component_mw = component_nodes[component_nodes['type']=='source']['capacity_mw'].sum()
        component_node_targets = component_nodes[component_nodes['type']=='target']
        total_component_gdp = component_node_targets['gdp'].sum()

        nominal_dict[ii] = {'nominal_mw': total_component_mw, 'nominal_gdp': total_component_gdp}

        if total_component_gdp != 0:
            component_target_mw_allocation = {target_id: total_component_mw*target_gdp/total_component_gdp for target_id, target_gdp in zip(component_node_targets.id.values, component_node_targets.gdp.values)}  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target has gdp:  tot_mw * target_gdp / tot_gdp. This is the nominal values
            nodes['nominal_mw'] = nodes['nominal_mw'] + nodes['id'].map(component_target_mw_allocation).fillna(0)  # maps the nominal values to the node dataframe

    if len(edges) != 0 and len(nodes) != 0:
        ## Split damaged components ##
        edges_damaged = edges[~edges.link.isin(edges_affected.link.values)]  # remove edges which are affected (edges_affected)
        G = create_graph(nodes, edges_damaged)  # new graph

        components = list(nx.connected_components(G))

        nodes['post_storm_mw'] = 0  # base level of mw value after storm
        for component in components:
            component_nodes = nodes[nodes['id'].isin(component)]
            total_component_mw_storm = component_nodes[component_nodes['type']=='source']['capacity_mw'].sum()
            component_node_targets = component_nodes[component_nodes['type']=='target']
            total_component_gdp_storm = component_node_targets['gdp'].sum()
            if total_component_gdp_storm != 0:

                ## Rerouting (method not yet verified) ##
                #component_target_mw_storm_allocation = {target_id: total_component_mw_storm*target_gdp/total_component_gdp_storm for target_id, target_gdp in zip(component_node_targets.id.values, component_node_targets.gdp.values)}  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target then has a new damaged mw: tot_mw_subnetwork * target_gdp / tot_gdp_subnetwork (for the sub-network in which the target is located)
                ## No Rerouting ##
                component_target_mw_storm_allocation = {target_id: total_component_mw_storm*target_gdp/nominal_dict[jj]['nominal_gdp'] for target_id, target_gdp, jj in zip(component_node_targets.id.values, component_node_targets.gdp.values, component_node_targets.component.values)}  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target then has a new damaged mw: tot_mw_subnetwork * target_gdp / tot_gdp (for the sub-network in which the target is located).

                nodes['post_storm_mw'] = nodes['post_storm_mw'] + nodes['id'].map(component_target_mw_storm_allocation).fillna(0)  # maps the nominal values to the node dataframe

        nodes['mw_loss_storm'] = nodes['nominal_mw'] - nodes['post_storm_mw']  # calculate mw loss
        nodes['f_value'] = 1 - nodes['mw_loss_storm'] / nodes['nominal_mw']  # calculate f value: power_after_storm / nominal_power
        nodes['gdp_damage'] = (1 - nodes['f_value']) * nodes['gdp']  # equivalent gdp value


else:
    print(f"No data in windfile for {nh}")
    nodes = pd.DataFrame()

print(f"{nh}: - saving")

storm_path = os.path.join(
    "data", "intersection", "storm_data", "individual_storms", region, sample, f"storm_{nh}"
)
if not os.path.exists(storm_path):
    os.makedirs(storm_path)

if len(edges_affected) != 0:  # to prevent writing empty dataframe
    edges_affected.to_file(
        os.path.join(
            storm_path, f"edges_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
        ),
        driver="GPKG",
    )
if len(polys_affected) != 0:  # to prevent writing empty dataframe
    polys_affected.to_file(
        os.path.join(
            storm_path, f"units_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
        ),
        driver="GPKG",
    )


## Targets ##

targets = gpd.GeoDataFrame()

if len(nodes) != 0:
    for box_id in nodes.box_id.unique():
        targets = targets.append(gpd.read_file(os.path.join(
                        "data",
                        "processed",
                        "all_boxes",
                        f"{box_id}",
                        f"targets_{box_id}.gpkg",)))  # add target of box

    targets = target_mapper('f_value', targets, nodes)  # map f_value
    targets = target_mapper('mw_loss_storm', targets, nodes)  # map mw loss after storm
    targets = target_mapper('gdp_damage', targets, nodes)  # map gdp damage from storm





target_cols = [
    "index",
    "area_km2",
    "population",
    "population_density_at_centroid",
    "gdp_pc",
    "gdp",
    "type",
    "box_id",
    "id",
    "f_value",
    "mw_loss_storm",
    "gdp_damage",
    "geometry",
]
if len(targets) == 0:  # if empty
    targets = gpd.GeoDataFrame(columns=target_cols)
    targets.loc[0, :] = [None] * len(target_cols)

targets.to_file(
    os.path.join(storm_path, f"targets__storm_r{region}_s{sample}_n{nh}.gpkg"),
    driver="GPKG",
)


# write storm track file
if len(TC) != 0:
    print(f"- writing {nh} to storm track file")
    TC_nh = TC[TC["number_hur"] == nh]

    #mask = TC["lon"] > 180.0
    #TC["lon"][mask] = TC["lon"] - 360.0  # adjust long

    TC_nh['lon'] = TC_nh['lon'].apply(lambda x: x if x <= 180 else x - 360)

    coords = [((lon, lat)) for lon, lat in zip(TC_nh["lon"], TC_nh["lat"])]
    storm_track = gpd.GeoDataFrame({"geometry": [LineString(coords)]})
    storm_track.to_file(
        os.path.join(storm_path, f"storm_track_r{region}_s{sample}_n{nh}.gpkg"),
        driver="GPKG",
    )



#%% add stats
today = date.today()


if not isNone(targets):
    f_75_1_temp, f_50_75, f_25_50, f_0_25_temp  = targets['f_value'].value_counts(bins=[0, 0.25, 0.5, 0.75, 1]).values.astype(float)  # note order
    f_0 = len(targets[targets['f_value']==0])
    f_0_25 = f_0_25_temp - f_0
    f_1 = len(targets[targets['f_value']==1])

    f_75_1 = f_75_1_temp - f_1
    print(f"f_75_1: {f_75_1}, f_75_1_temp: {f_75_1_temp}, f_1: {f_1}")

    totdamage = targets.gdp_damage.sum()

    num_affected = len(targets) - f_1
else:
    f_0, f_0_25, f_25_50, f_50_75, f_75_1 = [0]*5
    num_affected = 0
    totdamage = 0

stats_add = {
    "Storm ID": [nh],
    "Storm Region": [region],
    "Damages (gdp)": [totdamage],
    "targets affected": [num_affected],
    "targets operational 100%>op>75%": [f_75_1],
    "targets operational 75%>=op>50%": [f_50_75],
    "targets operational 50%>=op>25%": [f_25_50],
    "targets operational 25%>=op>0%": [f_0_25],
    "targets not operational (op=0%)": [f_0],
    "sim_run_date": [today.strftime("%d/%m/%Y")],
}

damagescsvpath = os.path.join(
    storm_path,
    f"storm_r{region}_s{sample}_n{nh}.txt",
)
with open(
    damagescsvpath, "w"
) as stormfile:  # open (overwrite) file for each storm year
    json.dump(stats_add, stormfile)
