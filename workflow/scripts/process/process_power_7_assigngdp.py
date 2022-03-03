"""Assigns each edge its gdp flow and notes the source-sink paths passing through it"""


# TODO
import sys

if 'linux' not in sys.platform:
    box_id = "box_1696"
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
else:
    box_id = sys.argv[1]

from process_power_functions import adj, adjbox, read_network, read_simple_network
from importing_modules import *

export_full_subgraph = False  # if True, will export an addional file for network which contains all adjacent connections (larger file) - Recommend False
include_paths = False  # if True will include paths running through edge (in the gpkg file - larger file) - Recommend False


def long2short(x):
    """Converts from long notation to short notation (more efficient to work with). x is string"""
    return x.replace('intermediate_','i').replace('_box_','b').replace('target_','t').replace('source_','s').replace('conn_','c')

def short2long(x):
    """Converts from short notation back to long notation. Order is crucial to avoid incorrect replacements. x is string"""
    return x.replace('c','conn_').replace('s','source_').replace('t','target_').replace('b','_box_').replace('i','intermediate_')

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
        # print(f"Examine {base_box} with boxes adjacent: {test_boxes}")
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
                # print(f"cant find box_{test_box} (nofile)")
                test_boxes.remove(test_box)  # remove since no more connections
                continue  # do not continue
            if len(connector_adj) == 0:
                # print(f"cant find box_{test_box} (len-0)")
                test_boxes.remove(test_box)  # remove since no more connections
                continue

            # print(f"loading box_{test_box}")
            test_fname = os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"box_{test_box}",
                f"simple_network_box_{test_box}.gpkg",
            )
            try:
                test_box_base_connections = connector_adj[f"box_{base_box}"]
            except:
                # print(f"No connections from {base_box} to {test_box}")
                test_boxes.remove(test_box)
                continue
            if len(test_box_base_connections) == 0:
                # print(f"No connections from {base_box} to {test_box} (len-0)")
                test_boxes.remove(test_box)
                continue

            try:
                test_nodes, test_edges = read_simple_network(test_fname)
            except:
                print(f"Simple network box does not exist for {test_box}. Will assume no connections between {base_box} and {test_box}.\nThis message should disappear with all boxes processed.")  # occurs when simple boxes not completed everywhere.
                test_boxes.remove(test_box)
                continue


            # print(f"Connections from {base_box} to {test_box}")
            # add adjacent box network

            test_G = create_graph_simple(test_nodes, test_edges)
            test_components = list(nx.connected_components(test_G))
            test_components_keep = []  # indices of test_components to include
            for (
                conn
            ) in test_box_base_connections:  # for each box connection to base_box
                if conn[5] not in list(nodes["id"]):
                    continue
                test_box_conn = conn[4]  # to_id
                for conn_comp in test_components:  # for each subgraph in test_box
                    if (
                        test_box_conn in conn_comp
                    ):  # if it contains a connection to base_box
                        # print(f"Found connections to base_box")
                        test_components_keep += conn_comp  # add it to the list

            test_nodes_keep = test_nodes[
                test_nodes["id"].isin(test_components_keep)
            ]  # keep only subgraphs connected to base_box
            test_edges_keep = test_edges[
                np.maximum(
                    test_edges["from_id"].isin(test_components_keep),
                    test_edges["to_id"].isin(test_components_keep),
                )
            ]  # keep only subgraphs connected to base_box either from_id or to_id (note max(True, False) = True)

            nodes = nodes.append(test_nodes_keep)  # add to network
            edges = edges.append(test_edges_keep)  # add to network

            test_box_base_connections_keep = [
                item
                for item in test_box_base_connections
                if item[5] in list(nodes["id"])
            ]  # only keep if in the subgraph

            if len(test_box_base_connections_keep) >= 1:
                # add adjacent box network connections
                all_connections_lst = [
                    [
                        connection_det[i]
                        for connection_det in test_box_base_connections_keep
                    ]
                    for i in range(len(test_box_base_connections_keep[0]))
                ]  # reorganise connection details
                all_connections_lst[-2] = [
                    sw.loads(str_ls) for str_ls in all_connections_lst[-2]
                ]  # convert to geom
                update_dict = dict(zip(edge_cols, all_connections_lst))
                edges = edges.append(gpd.GeoDataFrame(update_dict))  # add to network

                # print(f"added connections, adding {test_box} to to_check")
                to_test.append(
                    test_box
                )  # since test_box is connected, note to test all around test_box too

        # print(f"removing {base_box} from to_test\n")
        to_test.remove(base_box)  # tested all surrounding

    return nodes, edges, links_in_base_box


def create_graph(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from(
        (n.id, {"id": n.id, "type": n.type, "gdp": n.gdp, "capacity_mw": n.capacity_mw})
        for n in nodes.itertuples()
    )
    G.add_edges_from(
        (e.from_id, e.to_id, {"id": e.id, "length_m": e.length_eff})  # effective length is
        for e in edges.itertuples()
    )
    return G


def create_graph_simple(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from((n.id) for n in nodes.itertuples())
    G.add_edges_from((e.from_id, e.to_id) for e in edges.itertuples())
    return G


def assign_node_edge_gdp(G, links_in_base_box):
    components = list(nx.connected_components(G))
    node_gdps = []
    edge_gdps = []
    component_lst = []
    edge_gdp_sorted = {}

    for component in tqdm(components, desc="assigning gdp", total=len(components)):
        (
            c_node_gdp,
            c_edge_gdp,
            c_component,
            edge_gdp_indiv,
        ) = assign_component_gdp(
            G, component, links_in_base_box
        )  # c_components returns sources and sinks for each edge
        if len(c_node_gdp):
            node_gdps.append(c_node_gdp)
        if len(c_edge_gdp):
            edge_gdps.append(c_edge_gdp)
        if len(c_component):
            component_lst.append(c_component)

        edge_gdp_sorted = {
            **edge_gdp_sorted,
            **edge_gdp_indiv,
        }  # add new dictionary on the end. Each edge contains the target within the flows of itself

    def g0(df):
        if len(df) > 0:
            return pd.concat(df)
        else:
            return df

    node_gdp = g0(node_gdps)
    edge_gdp = g0(edge_gdps)
    component_df = g0(component_lst)

    return node_gdp, edge_gdp, component_df, edge_gdp_sorted


def edge_link_ids_from_nodes(G, route_nodes):
    next_nodes = iter(route_nodes)
    next(next_nodes)
    return ["__".join(sorted([u, v])) for u, v in zip(route_nodes, next_nodes)]

#@profile
def assign_component_gdp(G, component, links_in_base_box):
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
            return n.capacity_mw * c_gdp / c_cap
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


    print(f"{box_id} -- Sources: {len(sources)}, targets: {len(targets)}")
    if len(sources) and len(targets):
        for u in tqdm(sources.itertuples(), desc="sources", total=len(sources)):
            paths = nx.shortest_path(c, source=u.id, weight='length_eff')  # TODO add weight length_m
            for v in targets.itertuples():
                path = paths[v.id]  # path is the route from source to target
                path_gdp = u.gdp * (
                    v.gdp / c_gdp
                )  # (source gdp) * ( (target gdp)/(total component gdp) ) = (target gdp) * ( (source MW) / (total MW) )

                for link_id in edge_link_ids_from_nodes(c, path):

                    #if link_id in links_in_base_box:  # only save if in the base_box (to be examined)
                    edge_links[
                        link_id
                    ] += path_gdp  # each edge is given the associated gdp

                    # if include_paths:  # only do if include in gpkg
                    #     comp_path[link_id].append(
                    #         [path[0], path[-1]]
                    #     )  # [source, target]

                    route_id = (
                        path[0] + "_" + path[-1]
                    )  # source_sink unique for a flow
                    if link_id not in edge_gdp_indiv:
                        edge_gdp_indiv[
                            link_id
                        ] = (
                            {}
                        )  # add link_id dict to edge_gdp_indiv if not yet there

                    if route_id in edge_gdp_indiv[link_id]:
                        edge_gdp_indiv[link_id][
                            route_id
                        ] += path_gdp  # add path_gdp to target (source unknown) if target in edge_gdp_indiv[link_id]
                    else:
                        edge_gdp_indiv[link_id][
                            route_id
                        ] = path_gdp  # create path_gdp for target (source unknown) if target not in edge_gdp_indiv[link_id]

    print('Shortening base component names')
    links_in_base_box_shortened = [long2short(x) for x in links_in_base_box]
    print('filter base box...')
    edge_gdp_indiv = {key_keep: edge_gdp_indiv[key_keep] for key_keep in edge_gdp_indiv.keys() if key_keep in links_in_base_box_shortened}  # just keep in components in the base_box
    print('...filtered')
    #edge_gdp_indiv = {short2long(k): {short2long(kk): vv for kk,vv in v.items()} for k,v in edge_gdp_indiv.items()}  # change values back to original

    c_edges = pd.DataFrame({"link": k, "gdp": v} for k, v in edge_links.items())

    c_components = pd.DataFrame({"link": k, "path": v} for k, v in comp_path.items())

    return c_nodes[["id", "gdp"]], c_edges, c_components, edge_gdp_indiv


#%% run
s = time.time()
fname = os.path.join("data", "processed", "all_boxes", box_id, f"network_{box_id}.gpkg")

print("finding boxes within region")
print(f"{box_id}: importing and stitching together")
nodes, edges, links_in_base_box = combine_networks(box_id)
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
    comp_path,
    edge_gdp_sorted,
) = assign_node_edge_gdp(G, links_in_base_box)
# print(f"time: {round((time.time()-s)/60,2)}")


print(f"{box_id} -- saving edge_gdps_sorted")
with open(
    os.path.join(
        "data", "processed", "all_boxes", box_id, f"edge_gdp_sorted_{box_id}.txt"
    ),
    "w",
) as sortedjson:
    json.dump(edge_gdp_sorted, sortedjson)

out_fname = fname.replace("network", "network_with_gdp")


# print("node processing")
nodes = nodes.drop("gdp", axis=1)
if len(node_gdp) != 0:
    nodes = nodes.merge(node_gdp, on="id")
# print(f"time: {round((time.time()-s)/60,2)}")

# print("merging edges")
if len(comp_path) != 0:
    edge_gdp = edge_gdp.merge(
        comp_path, on="link"
    )  # merge the component sink source info for each edge
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
