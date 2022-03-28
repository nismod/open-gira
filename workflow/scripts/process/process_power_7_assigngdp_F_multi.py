"""Assigns each edge its gdp flow and notes the source-sink paths passing through it"""


# TODO
import sys

import pandas as pd
from pathos.multiprocessing import ProcessPool, cpu_count
from functools import partial
nodesuse = cpu_count() -2
if 'linux' not in sys.platform:
    box_id = "box_1794"  # 811, 810, 740, X 738, 668, 739, 666, 669, 667,
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
    nodesuse = 8
else:
    box_id = sys.argv[1]



from process_power_functions import adj, adjbox, read_network, read_simple_network, long2short, short2long
from importing_modules import *
import h5py
import pickle

export_full_subgraph = False  # if True, will export an addional file for network which contains all adjacent connections (larger file) - Recommend False
include_paths = False  # if True will include paths running through edge (in the gpkg file - larger file) - Recommend False



def flow_file(link_id):
    box_id = f"box_{link_id.split('b')[-1]}"
    return os.path.join('data', 'processed', 'tester', box_id+"_"+link_id+'.pkl')

#@profile
def combine_networks(
    box_id_orig,
):
    """finds adjacent boxes which are connected to original box_id and adds to nodes and edges.
    Returns all nodes and edges (both connected into base_box) and a list of links in the base_box"""
    fname = os.path.join(
        "data", "processed", "all_boxes", box_id_orig, f"network_{box_id_orig}.gpkg"
    )
    nodes, edges = read_network(fname)

    nodes["orig"] = True  # note that these are from the base box
    edges["orig"] = True  # note that these are from the base box

    edges['length_eff'] = edges.geometry.length  # add the geometry lengths

    links_in_base_box = list(edges["link"])  # notes the links in the base box
    if links_in_base_box[0] != None:
        links_in_base_box = set([long2short(x) for x in links_in_base_box])  # shorten names
    else:
        links_in_base_box = set()


    to_test = [int(box_id_orig[4:])]  # start with
    edge_cols = list(edges.columns)
    checked = []  # list of checked boxes (to note not to recheck)
    while len(to_test) >= 1:
        if len(checked) % 1  == 0:
            print(f"{box_id_orig}: checked {len(checked)} boxes")

        to_test = list(set(to_test))  # get rid of duplicates

        base_box = to_test[0]
        test_boxes = adj(base_box)  # boxes around base_box
        test_boxes_copy = test_boxes.copy()
        print(f"Examine {base_box} with boxes adjacent: {test_boxes}")
        for test_box in test_boxes_copy:
            if (
                f"{min(base_box, test_box)}_{max(base_box, test_box)}" in checked
            ):  # ensures no re-checks
                continue
            checked.append(f"{min(base_box, test_box)}_{max(base_box, test_box)}")
            try:  # try and open the network file
                with open(
                    os.path.join(
                        "data",
                        "processed",
                        "all_boxes",
                        f"box_{test_box}",
                        f"connector_box_{test_box}.txt",
                    ),
                    "r",
                ) as file_ex:
                    connector_adj = json.load(file_ex)
            except:  # if doesnt exist
                print(f"cant find box_{test_box} (nofile)")
                test_boxes.remove(test_box)  # remove since no more connections
                continue  # do not continue
            if len(connector_adj) == 0:
                print(f"cant find box_{test_box} (len-0)")
                test_boxes.remove(test_box)  # remove since no more connections
                continue

            print(f"loading box_{test_box}")
            test_fname = os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"box_{test_box}",
                f"network_box_{test_box}.gpkg",
            )
            try:
                test_box_base_connections = connector_adj[f"box_{base_box}"]
            except:
                print(f"No connections from {base_box} to {test_box}")
                test_boxes.remove(test_box)
                continue
            if len(test_box_base_connections) == 0:
                print(f"No connections from {base_box} to {test_box} (len-0)")
                test_boxes.remove(test_box)
                continue

            try:
                test_nodes, test_edges = read_network(test_fname)
            except:
                print(f"network box does not exist for {test_box}. Will assume no connections between {base_box} and {test_box}.\nThis message should disappear with all boxes processed.")  # occurs when simple boxes not completed everywhere.
                test_boxes.remove(test_box)
                continue


            print(f"Connections from {base_box} to {test_box}")
            # add adjacent box network

            # test_G = create_graph(test_nodes, test_edges)
            # test_components = list(nx.connected_components(test_G))
            # test_components_keep = []  # indices of test_components to include
            # for (
            #     conn
            # ) in test_box_base_connections:  # for each box connection to base_box
            #     if conn[5] not in list(nodes["id"]):  # if not already included
            #         continue
            #     test_box_conn = conn[4]  # to_id
            #
            #     #print(test_nodes[test_nodes['id']==test_box_conn]['component_id'], test_box)
            #     try:
            #         component_id = int(test_nodes[test_nodes['id']==test_box_conn]['component_id'])
            #         if component_id not in test_components_keep:  # if not already in list
            #             test_components_keep.append(component_id)  # add it to the list
            #     except:
            #         print(f'{base_box} excepted')
            #         pass

            test_G = create_graph(test_nodes, test_edges)
            test_components = list(nx.connected_components(test_G))
            test_components_keep = []  # indices of test_components to include
            for (
                conn
            ) in test_box_base_connections:  # for each box connection to base_box (from_id in base_box)
                print('here:' ,conn['to_id'])
                #if not conn['to_id']:  # if not empty
                if conn['to_id'] not in list(nodes["id"]):   # if adjacent box network connection is not in the node df, continue
                    continue
                test_box_conn = conn['from_id']  #
                for conn_comp in test_components:  # for each subgraph in test_box
                    if (
                        test_box_conn in conn_comp
                    ):  # if it contains a connection to base_box
                        # print(f"Found connections to base_box")
                        test_components_keep += conn_comp  # add it to the list




            # test_nodes_keep = test_nodes[
            #     test_nodes["component_id"].isin(test_components_keep)
            # ]  # keep only subgraphs connected to base_box

            test_nodes_keep = test_nodes[
                test_nodes["id"].isin(test_components_keep)]

            if len(test_edges) != 1 and test_edges['id'].iloc[0] != None:
                test_edges_keep = test_edges[
                    np.maximum(
                        test_edges["from_id"].isin(test_nodes_keep['id']),
                        test_edges["to_id"].isin(test_nodes_keep['id']),
                    )
                ]  # keep only subgraphs connected to base_box either from_id or to_id (note max(True, False) = True)

                #print(nodes['id'].iloc[0][-4:])
                # test = '0'
                # try:
                #     test = test_nodes_keep['id'].iloc[0][-4:]
                # except:
                #     a = 1
                # if test == '1650':
                #     print('---------------> hit')
                print('adding nodes/edges')
                nodes = nodes.append(test_nodes_keep)  # add to network  # TODO why counted double
                edges = edges.append(test_edges_keep)  # add to network

            test_box_base_connections_keep = [
                item
                for item in test_box_base_connections
                if item['to_id'] in list(nodes["id"])
            ]  # only keep if in the subgraph

            if len(test_box_base_connections_keep) >= 1:
                # add adjacent box network connections
                # all_connections_lst = [
                #     [
                #         connection_det[i]
                #         for connection_det in test_box_base_connections_keep
                #     ]
                #     for i in range(len(test_box_base_connections_keep[0]))
                # ]  # reorganise connection details
                # test_box_base_connections_keep['geome'] = [
                #     sw.loads(str_ls) for str_ls in all_connections_lst[-2]
                # ]  # convert to geom
                #update_dict = dict(zip(edge_cols, test_box_base_connections_keep))
                update_gdf = gpd.GeoDataFrame(test_box_base_connections_keep)
                update_gdf['geometry'] = update_gdf['geometry'].apply(sw.loads)
                edges = edges.append(update_gdf)  # add to network

                print(f"added connections, adding {test_box} to to_check")
                to_test.append(
                    test_box
                )  # since test_box is connected, note to test all around test_box too


            print(len(nodes))

        print(f"removing {base_box} from to_test\n")
        to_test.remove(base_box)  # tested all surrounding

    nodes = nodes.drop_duplicates()  # catch any errors  # TODO
    edges = edges.drop_duplicates()  # catch any errors   # TODO

    #nodes.to_file("testme.gpkg", layer="nodes", driver="GPKG")   # TODO
    #edges.to_file("testme.gpkg", layer="edges", driver="GPKG")  # TODO
    return nodes, edges


def create_graph(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from(
        (n.id, {"id": n.id, "type": n.type, "gdp": n.gdp, "capacity_mw": n.capacity_mw})
        for n in nodes.itertuples()
    )
    G.add_edges_from(
        (e.from_id, e.to_id, {"id": e.id, "length_m": e.geometry.length})  # effective length is
        for e in edges.itertuples()
    )
    return G



def create_graph_simpleXX(nodes, edges):
    G = nx.Graph()
    if len(nodes) != 1 and nodes['id'].iloc[0] != None:
        G.add_nodes_from((n.id) for n in nodes.itertuples())
    if len(edges) != 1 and edges['id'].iloc[0] != None:
        G.add_edges_from((e.from_id, e.to_id) for e in edges.itertuples())
    return G


def assign_node_edge_gdp(G):
    components = list(nx.connected_components(G))
    node_gdps = []
    edge_gdps = []
    component_lst = []
    route_gdp = {}

    for component in tqdm(components, desc="assigning gdp", total=len(components)):
        (
            c_node_gdp,
            c_edge_gdp,
            route_gdp_indiv,
            ) = assign_component_gdp(
            G, component
        )  # c_components returns sources and sinks for each edge
        if len(c_node_gdp):
            node_gdps.append(c_node_gdp)
        if len(c_edge_gdp):
            edge_gdps.append(c_edge_gdp)
        route_gdp = {**route_gdp, **route_gdp_indiv}


    def g0(df):
        if len(df) > 0:
            return pd.concat(df)
        else:
            return df

    node_gdp = g0(node_gdps)
    edge_gdp = g0(edge_gdps)

    return node_gdp, edge_gdp, route_gdp


def edge_link_ids_from_nodes(G, route_nodes):
    next_nodes = iter(route_nodes)
    next(next_nodes)
    return ["__".join(sorted([u, v])) for u, v in zip(route_nodes, next_nodes)]


def target_process(v, u, paths, c_gdp, c):  # TODO possinly directly use partial to inpout edge_routes etc dict to directly manuplaute

    edge_routes = dict()
    edge_links_gdp = defaultdict(int)
    path = paths[v.id]  # path is the route from source to target
    path_gdp = u.gdp * (
        v.gdp / c_gdp
    )  # (source gdp) * ( (target gdp)/(total component gdp) ) = (target gdp) * ( (source MW) / (total MW) )


    for link_id in edge_link_ids_from_nodes(c, path):

        edge_links_gdp[
            link_id
        ] += path_gdp  # each edge is given the associated gdp

        if link_id not in edge_routes.keys():
            edge_routes[link_id] = {(path[0], path[-1])}
        else:
            edge_routes[link_id].add((path[0], path[-1]))

    source_sink_gdp = [v.id, path_gdp]

    return edge_routes, source_sink_gdp, edge_links_gdp

#@profile
def assign_component_gdp(G, component):
    # create connected component subgraph
    c = G.subgraph(component).copy()
    c = nx.relabel_nodes(c, lambda x: long2short(x))  # relabel for efficiency
    # nodes dataframe
    c_nodes = pd.DataFrame(n for _, n in c.nodes(data=True))

    # component total gdp and capacity
    assert len(c_nodes.columns) > 0
    c_gdp = c_nodes.gdp.sum()
    c_cap = c_nodes.capacity_mw.sum()

    # assign GDP to nodes in proportion to capacity
    def assign_source_gdp(n):
        if n.type == "source":
            return n.capacity_mw * c_gdp / c_cap  # floats added for case when 0 is a string
            # return float(n.capacity_mw) * float(c_gdp) / float(c_cap)  # floats added for case when 0 is a string
        return n.gdp

    c_nodes["gdp"] = c_nodes.apply(assign_source_gdp, axis=1)  # TODO check midpoints not affected

    # assign GDP "flow" along shortest path from source to target, sharing source
    # gdp proportionally between targets
    edge_links = defaultdict(int)
    comp_path = defaultdict(list)  # keeps track of source and sink nodes for each edge

    sources = c_nodes[c_nodes.type == "source"].copy()
    targets = c_nodes[c_nodes.type == "target"].copy()

    sources['id'] = sources['id'].str.replace('source_','s').str.replace('_box_','b')  # improve computational efficiency
    targets['id'] = targets['id'].str.replace('target_','t').str.replace('_box_','b')  # improve computational efficiency

    edge_gdp_indiv = dict()

    edge_routes = dict()  # {edge: set_of_all_routes_through_that_edge}

    route_gdp_indiv_empty = dict(zip(sources['id'], [{}]*len(sources['id'])))  # {source1:{target1: gdp_flow1, target2: gdp_flow2, ...}, source2: ... }
    route_gdp_indiv = dict()
    #route_gdp_indiv_empty.copy()


    c_edges = pd.DataFrame(columns=['link', 'gdp'])


    base_flow_path = os.path.join('data', 'processed', 'all_boxes', box_id, 'gdp_flows')   # TODO for all affected
    if not os.path.exists(base_flow_path):
        os.makedirs(base_flow_path)


    print('starting pool')
    pool = ProcessPool(nodes=nodesuse)



    print(f"{box_id} -- Sources: {len(sources)}, targets: {len(targets)}")
    if len(sources) and len(targets):
        for u in tqdm(sources.itertuples(), desc="sources", total=len(sources)):

            paths = nx.shortest_path(c, source=u.id, weight='length_m')  # TODO add weight length_m


            print('starting partial')
            pool_partial=partial(target_process, u=u, paths=paths, c_gdp=c_gdp, c=c)#, edge_routes=edge_routes, route_gdp_indiv=route_gdp_indiv)
            print('starting output')
            output = pool.map(pool_partial, targets.itertuples())
            edge_routes_output = [x[0] for x in output]
            # include from output. Note that edge_routes IS reset after each source loop
            print('pool output complete')

            link_ids_u = set()  # set of all links affected by source u
            for edge_routes in edge_routes_output:
                link_ids_u.update(edge_routes.keys())  # add new link_ids
                for link_id, routes in edge_routes.items():
                    #if len(routes) != 0:
                    with open(flow_file(link_id), 'ab') as f:
                        _ = [pickle.dump(route, f) for route in routes]


            route_gdp_indiv_output = {x[1][0]:x[1][1] for x in output}  # {target1: gdp_u_to_target1, target2: gdp_u_to_target2, ...},
            route_gdp_indiv[u.id] = route_gdp_indiv_output  # include from output. Note that route_gdp_indiv is NOT reset after each source loop

            edge_links_gdp = [x[2] for x in output]


            c_edges = pd.merge(pd.DataFrame({'link':list(link_ids_u)}), c_edges, on='link', how='outer').fillna(0)  # , 'gdp':[0]*len(link_ids_u)

            for edge_link_gdp in edge_links_gdp:
                c_edges['gdp'] = c_edges['gdp']+c_edges['link'].map(edge_link_gdp) #TODO shouldnt be NANs


    return c_nodes[["id", "gdp"]], c_edges, route_gdp_indiv




# def midpoint_expand(source_target_dict, box_id):
#     """Takes source target dictionary with gdps and expands the midpoints according to the collapsed_sources_targets file"""
#     def collapse_file_name(box_id):
#         return os.path.join('data', 'processed', 'all_boxes', box_id, f'collapsed_sources_targets_{box_id}.txt')
#
#     for key, value in source_target_dict.items():
#         if 'midpoint' in key:
#     with open()




#%% run
if __name__ == "__main__":
    s = time.time()


      # TODO remove 2
    gdp_base_flows = os.path.join('data', 'processed', 'tester')
    if 'linux' in sys.platform:
        if os.path.exists(gdp_base_flows):
            print('removing files')
            files = glob.glob(os.path.join(gdp_base_flows, "*.pkl"))
            for file in files:
                os.remove(file)
        else:
            print('adding folder')
            os.makedirs(gdp_base_flows)
    else:  # windows
        if os.path.exists(gdp_base_flows):
            raise RuntimeError("Deleter gdp flows folder first")  # TODO
        os.makedirs(gdp_base_flows)





    fname = os.path.join("data", "processed", "all_boxes", box_id, f"network_{box_id}.gpkg")

    print("finding boxes within region")
    print(f"{box_id}: importing and stitching together")
    nodes, edges = combine_networks(box_id)
    print(f"{box_id}: finished stitching")
    # print(f"time: {round((time.time()-s)/60,2)}")

    if len(nodes) == 1 and nodes["id"].iloc[0] == None:
        nodes = gpd.GeoDataFrame(columns=nodes.columns)
    if len(edges) == 1 and edges["id"].iloc[0] == None:
        edges = gpd.GeoDataFrame(columns=edges.columns)


    # else:
    print(f"{box_id}: creating graph")
    G = create_graph(nodes, edges)
    # print(f"time: {round((time.time()-s)/60,2)}")


    # print("assigning node edges gdp")
    (
        node_gdp,
        edge_gdp,
        route_gdp,
    ) = assign_node_edge_gdp(G)
    # print(f"time: {round((time.time()-s)/60,2)}")

    a
    print(f"{box_id} -- saving edge_gdps_sorted")
    with open(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"edge_gdp_sorted_{box_id}.txt"
        ),
        "w",
    ) as sortedjson:
        json.dump(route_gdp, sortedjson)



    out_fname = fname.replace("network", "network_with_gdp")

    # print("node processing")
    nodes = nodes.drop("gdp", axis=1)
    if len(node_gdp) != 0:
        nodes = nodes.merge(node_gdp, on="id")
    # print(f"time: {round((time.time()-s)/60,2)}")


    if len(edge_gdp) != 0:
        #edges['link'] = [long2short(edge) for edge in edges['link']]
        edge_gdp['link'] = [short2long(edge) for edge in edge_gdp['link']]  # convert back to long form
        edges = edges.merge(edge_gdp, on="link")
        edges = edges[edges['box_id']==box_id]
    # print(f"time: {round((time.time()-s)/60,2)}")

    # print("cleaning data")
    edge_cols = list(edges.columns)
    if "path" in edge_cols:
        edges["path"] = [str(path) for path in edges["path"]]
    else:
        edge_cols.append("path")

    if len(edges) == 0:
        edges = gpd.GeoDataFrame(columns=edge_cols)
        edges.loc[0, :] = [None] * len(edge_cols)
    if len(nodes) == 0:
        nodes = gpd.GeoDataFrame(columns=nodes.columns)
        nodes.loc[0, :] = [None] * len(nodes.columns)

    if export_full_subgraph == True:
        print(f"{box_id}: Writing full subgraph files (in addition)")
        nodes.to_file(
            out_fname.replace("with_gdp", "with_gdp_full_subgraph"),
            layer="nodes",
            driver="GPKG",
        )
        edges.to_file(
            out_fname.replace("with_gdp", "with_gdp_full_subgraph"),
            layer="edges",
            driver="GPKG",
        )

    if len(edges) != 1 and edges["id"].iloc[0] != None:
        edges = edges[edges["orig"] == True]  # keep only in the base box
    if len(nodes) != 1 and nodes["id"].iloc[0] != None:
        nodes = nodes[nodes["orig"] == True]  # keep only in the base box
    edges.drop("orig", axis=1, inplace=True)
    nodes.drop("orig", axis=1, inplace=True)

    print(f"{box_id}: writing nodes to file")
    nodes.to_file(out_fname, layer="nodes", driver="GPKG")
    # print(f"time: {round((time.time()-s)/60,2)}")

    print(f"{box_id}: writing edges to file")
    edges.to_file(out_fname, layer="edges", driver="GPKG")

    print(f"\n{box_id} -- Time to run file: {(time.time() - s)/60:0.2f} minutes")
