"""Calculates the gdp losses for a given storm
"""

import os
import pandas as pd
import geopandas as gpd
import time
from datetime import date
from shapely.geometry import LineString, Point
import json
import sys
from damage_calculator import applythreshold
import networkx as nx
import numpy as np
import shapely.wkt as sw
from tqdm import tqdm
from geopy import distance


try:
    region = snakemake.params["region"]
    sample = snakemake.params["sample"]
    nh = snakemake.params["nh"]
    output_dir = snakemake.params['output_dir']
    reconstruction_cost_lowmedium = snakemake.params['reconstruction_cost_lowmedium']
    reconstruction_cost_high = snakemake.params['reconstruction_cost_high']
    central_threshold = snakemake.params['central_threshold']
    minimum_threshold = snakemake.params['minimum_threshold']
    maximum_threshold = snakemake.params['maximum_threshold']
    wind_file_start = snakemake.params['wind_file_start']
    wind_file_end = snakemake.params['wind_file_end']
    all_boxes = snakemake.params['all_boxes']
except:
    raise RuntimeError("Snakemake parameters not found")
    region = 'NA'
    sample = '0'
    nh = '0_148_11'
    output_dir = 'results'
    reconstruction_cost_lowmedium = 200000
    reconstruction_cost_high = 400000
    wind_file_start = 'STORM_DATA_CMCC-CM2-VHR4_'
    wind_file_end = '_IBTRACSDELTA'
    central_threshold = 43
    minimum_threshold = 39
    maximum_threshold = 47
    all_boxes = ['box_955', 'box_956', 'box_957', 'box_884']
    all_boxes = [f'box_{num}' for num in [809, 810, 811, 812, 880, 881, 882, 883, 884, 952, 955, 956, 957, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1097, 1098, 1099, 1103, 1104]]
    all_boxes = [f'box_{num}' for num in [884, 955, 956, 957, 1028, 1029, 1030, 1031, 1103, 1104]]





reconstruction_cost_high = float(reconstruction_cost_high)
assert reconstruction_cost_high >= 0
reconstruction_cost_lowmedium = float(reconstruction_cost_lowmedium)
assert reconstruction_cost_lowmedium >= 0
threshold_list = [central_threshold, minimum_threshold, maximum_threshold]

def isNone(df):
    """Checks if dataframe contains solely the None row (required for snakemake and gpkg files)"""
    if len(df) == 1 and df["id"].iloc[0] == None:
        return True
    else:
        return False


def read_edges_make_unique(fname):
    """Read edges and add unique link column"""
    edges = gpd.read_file(fname, layer="edges")

    if isNone(edges):  # no edges
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
    G.add_edges_from((e.from_id, e.to_id, {"id": e.link}) for e in edges.itertuples())
    return G


def box_connectors(box_id):
    """Loads and returns the connector dictionary"""
    with open(
        os.path.join(
            output_dir,
            "power_processed",
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
        output_dir, "power_processed", "all_boxes", box_id, f"network_{box_id}.gpkg"
    )
    return fname


def component_select(box_id, node_set, network_dict, all_boxes):
    """Finds the connected component(s) in box_id which contain(s) node_set
    Input:
        box_id: str
        node_set: set of nodes which are required to be in the connected component(s)
        network_dict: dictionary of (nodes, edges) corresponding to the box_id keys (to reduce time on loading, only do once)
        all_boxes: list of boxes outside of which not to be examined
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

    nodes_init, edges_init = network_dict[box_id]

    if isNone(nodes_init) or isNone(edges_init):  # if None row, return empty
        return nodes, edges, node_set, conn_boxes

    G_init = create_graph(nodes_init, edges_init)  # create a graph
    components_init = list(nx.connected_components(G_init))  # split into components
    for comp_init in components_init:  # for each component
        if any(
            x in comp_init for x in node_set
        ):  # if any node from the node_set and current node list is in this component

            nodes_comp = nodes_init[
                nodes_init["id"].isin(comp_init)
            ]  # nodes in this component
            edges_comp = edges_init[
                np.maximum(
                    edges_init["from_id"].isin(comp_init),
                    edges_init["to_id"].isin(comp_init),
                )
            ]  # edges in this component

            nodes = nodes.append(nodes_comp)  # keep the nodes of this component
            edges = edges.append(
                edges_comp
            )  # keep only edges of this component (note max(True, False) = True)

            ## note which nodes have been examined and which remain ##
            node_set_in_comp_init = {
                node for node in node_set if node in comp_init
            }  # nodes in component examined now, so will be removed
            conn_dict = {
                k: [
                    v_conn["to_id"]
                    for v_conn in v_lst
                    if v_conn["from_id"] in comp_init.difference(node_set_in_comp_init)
                ]
                for k, v_lst in connectors.items() if k in all_boxes
            }  # add the nodes (to_id) if the from_id is in the component (but not if connected to the previous box)
            node_set = node_set.difference(
                node_set_in_comp_init
            )  # node_set now contains all outer search nodes still to be examined
            conn_boxes = set(conn_dict.keys())  # note the boxes to still be examined
            conn_new = set().union(
                *conn_dict.values()
            )  # note the nodes still to be examined
            node_set = node_set.union(
                conn_new
            )  # update to include new node connections

            ## add the links for the examined cross-box connections ##
            cross_box_links = (
                []
            )  # list of dictionary with the relevant cross-box connections
            for connectors_indiv in connectors.values():
                cross_box_links += [
                    conn_link_dict
                    for conn_link_dict in connectors_indiv
                    if conn_link_dict["from_id"] in node_set_in_comp_init
                ]  # include the links where the 'from_id' is in the nodes in comp_init. That means the link travels into the box currently being examined
            if len(cross_box_links) != 0:
                update_gdf = gpd.GeoDataFrame(cross_box_links)
                update_gdf["geometry"] = update_gdf["geometry"].apply(
                    sw.loads
                )  # make geometry type
                edges = edges.append(update_gdf)  # add the new links

    return nodes, edges, node_set, conn_boxes


def combine_networks(
    edge_damaged, all_boxes
):
    """finds adjacent boxes which are connected to original box_id and adds to nodes and edges.
    Returns all nodes and edges (both connected into base_box). edge_damaged is dataframe input containing link and box_id and id
    Searches strictly only in all_boxes"""

    box_id_orig = edge_damaged.box_id
    fname = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id_orig, f"network_{box_id_orig}.gpkg"
    )
    nodes_orig, edges_orig = read_network(
        fname
    )  # load the network for the box which contains the damaged component
    G_orig = create_graph(nodes_orig, edges_orig)  # create a graph
    components_orig = list(nx.connected_components(G_orig))  # split into components
    for comp_orig in components_orig:  # for each component
        if any(
            x in comp_orig for x in [edge_damaged.from_id, edge_damaged.to_id]
        ):  # if either end (from_id or to_id) of the damaged component is in this component

            nodes = nodes_orig[
                nodes_orig["id"].isin(comp_orig)
            ]  # keep only the nodes of this component
            edges = edges_orig[
                np.maximum(
                    edges_orig["from_id"].isin(comp_orig),
                    edges_orig["to_id"].isin(comp_orig),
                )
            ]  # keep only edges of this component (note max(True, False) = True)

            break  # edge_damaged.id will not be in remaining components, stop looking

    ## For the original box, note the connectors to conn_set ##
    conn_dict = {
        key_box: {
            conn_dict["to_id"]
            for conn_dict in conn_lst
            if conn_dict["from_id"] in nodes.id.values
        }
        for key_box, conn_lst in box_connectors(box_id_orig).items()
    }  # if any link starts in the current nodes df, then add it to the dictionary

    conn_dict = {k: v for k, v in conn_dict.items() if k in all_boxes}  # filter


    conn_set = set().union(
        *conn_dict.values()
    )  # conn set will contain a set of links which have to be explored

    to_examine = set(conn_dict.keys())  # set of boxes to examine


    count = 0
    network_dict = (
        dict()
    )  # dictionary {box_id1: (nodes_of_box_id1, edges_od_box_id1), box_id2, ...}
    while len(conn_set) >= 1:  # while there are still links to be examined

        box_id_examine = min(
            to_examine
        )  # pick one from to_examine (which is irrelevant)

        #print(f'examining {box_id_examine}')
        if (
            box_id_examine not in network_dict.keys()
        ):  # if not in the dictionary, then add it
            network_dict[box_id_examine] = read_network(network_name(box_id_examine))

        count += 1
        if count % 100 == 0:
            print(f"Examined {count} boxes")
        to_examine = to_examine.difference({box_id_examine})

        nodes_examine, edges_examine, conn_set, to_examine_newboxes = component_select(
            box_id_examine, conn_set, network_dict, all_boxes
        )  # extract the network components which connect to any node in conn_set, then update conn_set to not include these are more
        # if len(nodes_examine) > 0:
        #     nodes_examine.to_file(os.path.join('tester', f'{count}_nodes.gpkg'), driver='GPKG')
        # if len(edges_examine) > 0:
        #     edges_examine.to_file(os.path.join('tester', f'{edge_damaged.id}_{count}_edges.gpkg'), driver='GPKG')

        nodes = nodes.append(nodes_examine)  # add
        edges = edges.append(edges_examine)  # add
        to_examine = to_examine.union(to_examine_newboxes)  # update

    return nodes, edges


def target_mapper(feature, targets, nodes):
    """Maps feature (of nodes df) to targets df. Will remove rows in which map has no value"""
    target_mapper_dict = {
        target_id: feature
        for target_id, feature, node_type in zip(
            nodes["id"].values, nodes[feature].values, nodes["type"].values
        )
        if node_type == "target"
    }  # dictionary {target1: f_value_of_target1, target2: f_value_of_target2, ... }
    targets[feature] = (
        targets["id"].map(target_mapper_dict).fillna("remove_me")
    )  # map to targets, note which are not connected (with 'remove_me')
    return targets[targets.f_value != "remove_me"]  # remove unwanted targets

def add_stats(targets, storm_path_thrval, direct_damage_cost):
    """Write the stats to a txt file. Requires target gpd, storm_path_thrval and direct damage cost value"""
    today = date.today()
    if not isNone(targets) and len(targets) != 0:
        f_0_25_temp, f_25_50, f_50_75, f_75_1_temp = (
            targets["f_value"]
            .value_counts(bins=[0, 0.25, 0.5, 0.75, 1], sort=False)
            .values.astype(float)
        )  # note order
        f_0 = len(targets[targets["f_value"] == 0])
        f_0_25 = f_0_25_temp - f_0
        f_1 = len(targets[targets["f_value"] == 1])

        f_75_1 = f_75_1_temp - f_1

        totdamage = targets.gdp_damage.sum()

        num_affected = len(targets) - f_1

        assert f_0 >= 0
        assert f_0_25 >= 0
        assert f_25_50 >= 0
        assert f_50_75 >= 0
        assert f_75_1 >= 0
        assert f_1 >= 0

        pop_affected = targets[targets['f_value']<1].population.sum()  # sum the population where the power after the storm is NOT the same as the nominal power (f<1 ie f!=1)
        pop_f0 = targets[targets['f_value']==0].population.sum()  # sum the population which has no power (f=0)
        pop_effective = ((1-targets['f_value'])*targets['population']).sum()  # effective population affected is the fraction of the population which has an effective power of 0 i.e. (1-f)*pop. This means with 100 people and f=0.2, pop_eff = 80.

        countries_affected = '_'.join(targets.country.unique())  # join to one simple string country1_country2_country3...  This greatly simplified json and pandas handling later

    else:
        f_0, f_0_25, f_25_50, f_50_75, f_75_1 = [0] * 5
        num_affected = 0
        totdamage = 0
        pop_affected = 0
        pop_f0 = 0
        pop_effective = 0
        countries_affected = None

    stats_add = {
        "Storm ID": [nh],
        "Storm Region": [region],
        "GDP losses": [totdamage],
        "targets affected": [num_affected],
        "targets 1>f>0_75": [f_75_1],
        "targets 0_75>=f>0_5": [f_50_75],
        "targets 0_5>=f>0_25": [f_25_50],
        "targets 0_25>=f>0": [f_0_25],
        "targets with no power (f=0)": [f_0],
        "population affected": [pop_affected],
        "population with no power (f=0)": [pop_f0],
        "effective population affected": [pop_effective],
        "affected countries": [countries_affected],
        "reconstruction cost": [direct_damage_cost],
        "sim_run_date": [today.strftime("%d/%m/%Y")],
    }

    damagescsvpath = os.path.join(storm_path_thrval, f"storm_r{region}_s{sample}_n{nh}.txt",
    )
    with open(
        damagescsvpath, "w"
    ) as stormfile:  # open (overwrite) file for each storm year
        json.dump(stats_add, stormfile)





def eval_coords(coords, type, wind, fragility_data):
    """Evaluate the coordinates and returns direct damage"""
    dist = 0
    fragility_curve = fragility_data[type]
    cost_per_km = fragility_curve['factor']*np.interp(wind, fragility_curve['curve']['wind_speed'], fragility_curve['curve']['fragility'])
    line_coords = [(x[1], x[0]) for x in coords]
    if len(line_coords) >= 2:
        for jj in range(len(line_coords)-1):
            dist += distance.distance(line_coords[jj], line_coords[jj+1]).km  # add the km length of that row
    cost = dist*cost_per_km

    return cost


def direct_damage(linestring_df):
    """For a dataframe where the geometry consists of linestrings, returns a list of direct damages (order of linestring_df) based on fragility curve"""


    damage_lst = []

    for ii in range(len(linestring_df)):
        damage = 0
        transmission_type = linestring_df.iloc[ii].source
        if transmission_type == None:
            transmission_type = 'openstreetmap'
        transmission_max_wind = linestring_df.iloc[ii].wind_location

        if type(linestring_df.iloc[ii].geometry) == type(LineString([[1,2],[3,4]])):  # check is linestring
            line_coords = list(linestring_df.iloc[ii].geometry.coords)  # extract the coordinates of a row
            damage += eval_coords(line_coords, transmission_type, transmission_max_wind, fragility_data)

        else:  #multistring
            for ms in range(len(linestring_df.iloc[ii]['geometry'].geoms)):
                line_coords = list(linestring_df.iloc[ii].geometry[ms].coords)  # extract the coordinates of a row

                damage += eval_coords(line_coords, transmission_type, transmission_max_wind, fragility_data)
        damage_lst.append(damage)

    return damage_lst


def map_col(df1, df2, col_name):
    '''Maps df on link from col_name on df1 to df2'''

    dict_map = dict(zip(df1['link'], df1[col_name]))
    df2[col_name] = df2['link'].map(dict_map)
    return df2


## End function ##
# open fragility curve now (once, rather than every time needed)
fragility_folder = os.path.join('workflow', 'scripts')
lowmedium_fragility = pd.read_csv(os.path.join(fragility_folder, 'lowmedium_fragility.csv'))
lowmedium_fragility.columns = ['wind_speed', 'fragility']
high_fragility = pd.read_csv(os.path.join(fragility_folder, 'high_fragility.csv'))
high_fragility.columns = ['wind_speed', 'fragility']
fragility_data = {'gridfinder': {'curve': lowmedium_fragility, 'factor':reconstruction_cost_lowmedium}, 'openstreetmap': {'curve': high_fragility, 'factor':reconstruction_cost_high}}  # dictionary, if type of key then use fragility curve of value








storm_path_base = os.path.join(
    output_dir,
    "power_intersection",
    "storm_data",
    "individual_storms",
    region,
    sample
)

storm_path = os.path.join(
    storm_path_base,
    f"storm_{nh}",
)
if not os.path.exists(storm_path):
    os.makedirs(storm_path)


# set iteration variables
routeid_damaged = (
    set()
)  # marks which source_sink routes already damaged -> do not double count. Form: {(source1, target1), (source2, target2), ...}
totdamage = 0  # total damage (ensuring no double counts)
targetsdamaged = {}  # if the operation value of the target is desired
edges_affected = gpd.GeoDataFrame()
targets = gpd.GeoDataFrame()
polys_affected = gpd.GeoDataFrame()


print(f"{nh}: loading data")
# print('loading tracks')
stormfile = os.path.join(
    output_dir,
    "input",
    "stormtracks",
    "events",
    f"{wind_file_start}{region}_1000_YEARS_{sample}{wind_file_end}.txt",
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
TC = TC[["year", "number", "lat", "lon", "radius", "cat", "wind"]]
TC["number_hur"] = (
    str(sample)
    + "_"
    + TC["year"].astype(int).astype(str)
    + "_"
    + TC["number"].astype(int).astype(str)
)

# print("Loading wind")
windfile = os.path.join(
    output_dir,
    "power_intersection",
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
    winds_ev_filtered = applythreshold(winds_ev_all, minimum_threshold)

    winds_ev_filtered_todict = winds_ev_filtered[['ID_point', 'wind_location']]
    winds_ev_filtered_todict = winds_ev_filtered_todict.groupby('ID_point').max().reset_index()  # pick max wind for each unit

    ID_affected = list(winds_ev_filtered_todict["ID_point"])
    ID2wind = dict(zip(winds_ev_filtered_todict['ID_point'], winds_ev_filtered_todict['wind_location']))  # dictionary {ID_point1: wind_speed1, ID_point2: wind_speed2, ...}
    box_id_affected = winds_ev_filtered["box_id"].unique()
    box_id_affected = [box for box in box_id_affected if box in all_boxes]  # only include if box is in examining boxes

    # print("- grid")
    grid_data = gpd.read_file(
        os.path.join(output_dir, "power_intersection", "regions", f"{region}_unit.gpkg")
    )

    polys_affected = grid_data[grid_data["ID_point"].isin(ID_affected)].rename(
        columns={"box_id": "box_id_poly"}
    )

    start = time.time()
    for jj, box_id in enumerate(
        box_id_affected
    ):  # extract the damaged edges using this for loop
        print(f"-- Inspecting for damage {jj+1}/{len(box_id_affected)} -- {box_id}")

        box_edges = gpd.read_file(
            os.path.join(
                output_dir, "power_processed",  "all_boxes", box_id, f"network_{box_id}.gpkg"
            ),
            layer="edges",
        )





        ### OPTION 1 - Just keep overlay on damaged ID points ###
        box_edges["link"] = box_edges.apply(
            lambda e: "__".join(sorted([e.from_id, e.to_id])), axis=1
        )  # consistent naming

        box_edges_affected_forid = box_edges.overlay(
            polys_affected, how="intersection"
        )  # keeps edges that are affected grid points (only a part has to be in)

        box_edges_affected = box_edges[box_edges.link.isin(box_edges_affected_forid['link'].unique())]



        # ### OPTION 2 - Any line that intersects an ID Point ###
        # box_edges["link"] = box_edges.apply(
        #     lambda e: "__".join(sorted([e.from_id, e.to_id])), axis=1
        # )  # consistent naming
        #
        # box_edges_affected_forid = box_edges.overlay(
        #     polys_affected, how="intersection"
        # )  # keeps edges that are affected grid points (only a part has to be in)
        #
        #
        # # print(f'len box_edges = {len(box_edges)}')
        # box_edges_affected = box_edges[box_edges.link.isin(box_edges_affected_forid['link'].unique())]
        # box_edges_affected = map_col(box_edges_affected_forid, box_edges_affected, 'ID_point')
        # box_edges_affected = map_col(box_edges_affected_forid, box_edges_affected, 'region')
        # box_edges_affected = map_col(box_edges_affected_forid, box_edges_affected, 'box_id_poly')




        ### END OPTION ###


        edges_affected = edges_affected.append(
            box_edges_affected
        )  # add to master list of damaged edges

    edges = gpd.GeoDataFrame(columns=["link"])
    nodes = gpd.GeoDataFrame()

    print("Starting network connection expansion")
    startx = time.time()

    for edge_damaged in tqdm(edges_affected.itertuples(), desc='iterating damaged edges', total=len(edges_affected)):

        if edge_damaged.link not in edges.link.values:
            # print('New link')
            # print('Searching Network')
            nodes_new, edges_new = combine_networks(edge_damaged, all_boxes)
            s1 = time.time()
            nodes = nodes.append(nodes_new)
            s2 = time.time()
            # print(f"Time for s1 is {s2 - s1}")
            edges = edges.append(edges_new)
            # print(f"Time for s2 is {time.time() - s2}")

    print(f"Search took {round((time.time() - startx)/60,1)} mins")

    ## Now all the damaged edges can be found in nodes & edges
    G = create_graph(nodes, edges)

    # if len(edges) > 0:
    #     edges.to_file(os.path.join('tester', f'ALL_edges.gpkg'), driver='GPKG')

    components = list(nx.connected_components(G))

    ## Set nominal values ##
    nodes["nominal_mw"] = 0  # set base to zero
    nominal_dict = (
        dict()
    )  # dictionary: {component1: {'nominal_mw': nominal_mw_1, 'nominal_gdp': nominal_gdp_1}, ... }
    for ii, component in enumerate(components):
        component_nodes = nodes[nodes["id"].isin(component)]
        nodes.loc[nodes.id.isin(component), "component"] = ii
        total_component_mw = component_nodes[component_nodes["type"] == "source"][
            "capacity_mw"
        ].sum()
        component_node_targets = component_nodes[component_nodes["type"] == "target"]
        total_component_gdp = component_node_targets["gdp"].sum()

        nominal_dict[ii] = {
            "nominal_mw": total_component_mw,
            "nominal_gdp": total_component_gdp,
        }

        if total_component_gdp != 0:
            component_target_mw_allocation = {
                target_id: total_component_mw * target_gdp / total_component_gdp
                for target_id, target_gdp in zip(
                    component_node_targets.id.values, component_node_targets.gdp.values
                )
            }  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target has gdp:  tot_mw * target_gdp / tot_gdp. This is the nominal values
            nodes["nominal_mw"] = nodes["nominal_mw"] + nodes["id"].map(
                component_target_mw_allocation
            ).fillna(
                0
            )  # maps the nominal values to the node dataframe


    #############################

    if len(edges_affected) != 0:
        edges_affected['wind_location'] = edges_affected['ID_point'].map(ID2wind)  # map the wind values for later filtering
        polys_affected['wind_location'] = polys_affected['ID_point'].map(ID2wind)



    for thrval in threshold_list:
        # redefinitions, other variables will be overwritten

        storm_path_thrval = os.path.join(storm_path, str(thrval))
        if not os.path.exists(storm_path_thrval):
            os.makedirs(storm_path_thrval)

        if 'wind_location' not in edges_affected.columns:
            print('No edges affected')
            add_stats(gpd.GeoDataFrame(), storm_path_thrval, 0)  # write stats to file
        else:


            edges_affected_thrval = edges_affected[edges_affected['wind_location']>=thrval]  # need to be redefined
            polys_affected_thrval = polys_affected[polys_affected['wind_location']>=thrval]  # need to be redefined


            if len(edges) != 0 and len(nodes) != 0:
                ## Split damaged components ##
                edges_damaged = edges[
                    ~edges.link.isin(edges_affected_thrval.link.values)
                ]  # remove edges which are affected (edges_affected)
                G = create_graph(nodes, edges_damaged)  # new graph

                components = list(nx.connected_components(G))

                nodes["post_storm_mw"] = 0  # base level of mw value after storm
                for component in components:
                    component_nodes = nodes[nodes["id"].isin(component)]
                    total_component_mw_storm = component_nodes[
                        component_nodes["type"] == "source"
                    ]["capacity_mw"].sum()
                    component_node_targets = component_nodes[
                        component_nodes["type"] == "target"
                    ]
                    total_component_gdp_storm = component_node_targets["gdp"].sum()
                    if total_component_gdp_storm != 0:

                        ## Rerouting (method not yet verified) ##
                        # component_target_mw_storm_allocation = {target_id: total_component_mw_storm*target_gdp/total_component_gdp_storm for target_id, target_gdp in zip(component_node_targets.id.values, component_node_targets.gdp.values)}  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target then has a new damaged mw: tot_mw_subnetwork * target_gdp / tot_gdp_subnetwork (for the sub-network in which the target is located)
                        ## No Rerouting ##
                        component_target_mw_storm_allocation = {
                            target_id: total_component_mw_storm
                            * target_gdp
                            / nominal_dict[jj]["nominal_gdp"]
                            for target_id, target_gdp, jj in zip(
                                component_node_targets.id.values,
                                component_node_targets.gdp.values,
                                component_node_targets.component.values,
                            )
                        }  # dictionary {target1: mw_for_target1, target2: mw_for_target2, ... }  Each target then has a new damaged mw: tot_mw_subnetwork * target_gdp / tot_gdp (for the sub-network in which the target is located).

                        nodes["post_storm_mw"] = nodes["post_storm_mw"] + nodes["id"].map(
                            component_target_mw_storm_allocation
                        ).fillna(
                            0
                        )  # maps the nominal values to the node dataframe

                nodes["mw_loss_storm"] = (
                    nodes["nominal_mw"] - nodes["post_storm_mw"]
                )  # calculate mw loss
                nodes["f_value"] = (
                    1 - nodes["mw_loss_storm"] / nodes["nominal_mw"]
                )  # calculate f value: power_after_storm / nominal_power
                nodes["gdp_damage"] = (1 - nodes["f_value"]) * nodes[
                    "gdp"
                ]  # equivalent gdp value




            print(f"{nh}, {thrval}m/s: - saving")


            #############################
            direct_damage_cost = 0
            if len(edges_affected_thrval) != 0:  # to prevent writing empty dataframe
                # calculate direct damage cost
                direct_damages = direct_damage(edges_affected_thrval)
                direct_damage_cost = sum(direct_damages)
                edges_affected_thrval['reconstruction_cost'] = direct_damages


                edges_affected_thrval.to_file(
                    os.path.join(
                        storm_path_thrval, f"edges_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
                    ),
                    driver="GPKG",
                )
            if len(polys_affected_thrval) != 0:  # to prevent writing empty dataframe
                polys_affected_thrval.to_file(
                    os.path.join(
                        storm_path_thrval, f"units_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
                    ),
                    driver="GPKG",
                )


            ## Targets ##

            targets = gpd.GeoDataFrame()
            if len(nodes) != 0:
                for box_id in nodes.box_id.unique():
                    targets = targets.append(
                        gpd.read_file(
                            os.path.join(
                                output_dir,
                                "power_processed",
                                "all_boxes",
                                f"{box_id}",
                                f"targets_{box_id}.gpkg",
                            )
                        )
                    )  # add target of box


                if 'target_3_box_957' in targets.id.values:
                    targets = targets[targets.id!='target_3_box_957']

                targets = target_mapper("f_value", targets, nodes)  # map f_value
                targets = target_mapper("mw_loss_storm", targets, nodes)  # map mw loss after storm
                targets = target_mapper("gdp_damage", targets, nodes)  # map gdp damage from storm


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
            if len(targets) != 0:  # if not empty
                targets.to_file(
                    os.path.join(storm_path_thrval, f"targets__storm_r{region}_s{sample}_n{nh}.gpkg"),
                    driver="GPKG",
                )



            add_stats(targets, storm_path_thrval, direct_damage_cost)  # write stats to file

    #############################
else:
    print(f"No data in windfile for {nh}")
    nodes = pd.DataFrame()



# write storm track file
if len(TC) != 0:
    print(f"- writing {nh} to storm track file")
    TC_nh = TC[TC["number_hur"] == nh].copy()

    TC_nh["lon"] = TC_nh["lon"].apply(lambda x: x if x <= 180 else x - 360)

    coords = [((lon, lat)) for lon, lat in zip(TC_nh["lon"], TC_nh["lat"])]

    if len(coords) >= 2:
        storm_track = gpd.GeoDataFrame({"geometry": [LineString(coords)]})  # stormtrack as a line
        storm_track.to_file(
            os.path.join(storm_path, f"storm_track_r{region}_s{sample}_n{nh}.gpkg"),
            driver="GPKG",
        )

        print('stormsize')
        storm_size = gpd.GeoDataFrame({"geometry":[Point(coord) for coord in coords], "radius":list(TC_nh.radius), "cat":list(TC_nh.cat), "wind":list(TC_nh.wind)})  # storm track points with extra features, for plotting
        storm_size.to_file(
            os.path.join(storm_path, f"storm_track_points_radius_r{region}_s{sample}_n{nh}.gpkg"),
            driver="GPKG",
        )
    else:
        print('No storm track data')
print(f'Finished_{nh}')
