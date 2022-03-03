"""Takes a complicated box with a network with connections in and out and simplifies it to just the in and out connections.

This file will take a box_id and generate the simplified network. If one does not wish to examine a box_id, but rather needs it (e.g. for assigning the gdp in a different(!) box_id), then one can use the simplified box. It is essentially a stand-in which drastically improves computational efficiency.

"""
import sys

if "linux" not in sys.platform:
    # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)


from process_power_functions import read_network
from importing_modules import *
import itertools as it

try:
    box_id = snakemake.params["box_id"]
    print(f"Running {box_id} box simplification")
except:
    box_id = "box_1942"  # TODO
    print('RUNNING FROM WINDOWS')


box_network_path = os.path.join(
    "data", "processed", "all_boxes", f"{box_id}", f"network_{box_id}.gpkg"
)
box_network_nodes, box_network_edges = read_network(box_network_path)
box_connector_path = os.path.join(
    "data", "processed", "all_boxes", f"{box_id}", f"connector_{box_id}.txt"
)

with open(box_connector_path, "r") as src:
    box_connector = json.load(src)
conns = []
border_nodes = []
for adj_box_dict in box_connector.values():
    conns += [x[6] for x in adj_box_dict if x != []]
    border_nodes += [x[4] for x in adj_box_dict if x != []]

x = [float(c.split(" ")[1][1:]) for c in conns]
y = [float(c.split(" ")[2][:-1]) for c in conns]
assert len(x) == len(y)
coords = [Point(x[i], y[i]) for i in range(len(x))]
border_nodes_cords = dict(zip(border_nodes, coords))


simple_file = os.path.join('data', 'processed', 'all_boxes', box_id, f'simple_network_{box_id}.gpkg')

if len(box_network_nodes) == 1 and box_network_nodes['id'].iloc[0] == None or len(box_network_edges) == 1 and box_network_edges['id'].iloc[0] == None:  # No connectors
    box_network_edges.to_file(simple_file, layer='edges', driver='GPKG')
    box_network_nodes.to_file(simple_file, layer='nodes', driver='GPKG')
else:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


        G = nx.Graph()
        G.add_nodes_from((n.id, {"type": n.type, "gdp": n.gdp, "capacity_mw": n.capacity_mw}) for n in box_network_nodes.itertuples())
        G.add_edges_from((e.from_id, e.to_id, {"length_m": e.geometry.length}) for e in box_network_edges.itertuples())
        components = list(nx.connected_components(G))

        component_midpoints = {}  # dictionary for midpoints of components
        box_simple_network_edges = gpd.GeoDataFrame()  # new data frame for simplified edges
        for ii, component in tqdm(enumerate(components), total=len(components), desc=f'processing components for {box_id}'):
            G_component = G.subgraph(component)  # create subgraph
            G_nodes = pd.DataFrame(n for _, n in G_component.nodes(data=True))
            tot_component_mw = G_nodes.capacity_mw.sum()  # total mw within that subgraph
            tot_component_gdp = G_nodes.gdp.sum()  # total gdp within that subgraph


            current_border_nodes = [
                border_node for border_node in border_nodes if border_node in component
            ]  # extract border nodes in this component

            if len(current_border_nodes) >= 1:
                if len(current_border_nodes) == 1:  # only one in/out
                    xloc = border_nodes_cords[current_border_nodes[0]].x
                    yloc = border_nodes_cords[current_border_nodes[0]].y
                    component_midpoints[ii] = {'geometry':Point(xloc, yloc), 'total_mw': tot_component_mw, 'total_gdp': tot_component_gdp}
                else:  # >= 2
                    node_combinations = list(
                        it.combinations(current_border_nodes, 2)
                    )  # list of tuples containing all connections within this component

                    all_x_sum = [
                        border_nodes_cords[node_combination[0]].x
                        + border_nodes_cords[node_combination[1]].x
                        for node_combination in node_combinations
                    ]  # sum x values for from and to nodes
                    all_y_sum = [
                        border_nodes_cords[node_combination[0]].y
                        + border_nodes_cords[node_combination[1]].y
                        for node_combination in node_combinations
                    ]  # sum x values for from and to nodes
                    component_midpoints[ii] = {'geometry': Point(
                        0.5 * sum(all_x_sum) / len(all_x_sum), 0.5* sum(all_y_sum) / len(all_y_sum)
                    ),
                    'total_mw': tot_component_mw,
                    'total_gdp': tot_component_gdp}
                    # add 'centroid' of all points (0.5 required because of summing both from and to node x values)


                    shortest_path_calculated = []  # list of presaved shortest path calculations
                    node_combinations.sort()


                    for node_combination in tqdm(node_combinations, desc=f'component combinations {box_id} {ii}', total=len(node_combinations)):
                        from_id = node_combination[0]
                        to_id = node_combination[1]

                        #shortest_length = nx.shortest_path_length(G_component, from_id, to_id, 'length_m')

                        if from_id not in shortest_path_calculated:  # if not already calculated then work out the shortest routes from the from_id node
                            shortest_length_all = nx.single_source_dijkstra_path_length(G_component, from_id, weight='length_m')  #  weight='length_m' <-- add this parameter for weight distance
                            shortest_path_calculated.append(from_id)

                        shortest_length = shortest_length_all[to_id]


                        new_line = {
                            "id": from_id + "__" + to_id,
                            "from_id": from_id,
                            "to_id": to_id,
                            "box_id": box_id,
                            "type": 'effective transmission line',
                            "geometry": LineString([border_nodes_cords[from_id], component_midpoints[ii]['geometry'], border_nodes_cords[to_id]]),  # line from from_id to to_id via midpoint
                            "length_eff": shortest_length
                        }  # add details
                        box_simple_network_edges = box_simple_network_edges.append(new_line, ignore_index=True)

        box_simple_network_nodes = pd.merge(box_network_nodes, pd.DataFrame({'id':border_nodes}), how='inner')  # keep only the nodes at the border

        for idx, midpoint in component_midpoints.items():  # for each midpoint add effective values which can be used to substitute the complex original network
            # midpoint
            midpoint_name = f'midpoint_{idx}_{box_id}'
            box_simple_network_nodes = box_simple_network_nodes.append({'id':midpoint_name, 'type':'midpoint', 'box_id':box_id, 'geometry':midpoint['geometry']}, ignore_index=True)  # add base midpoint

            # effective powerplant
            midpoint_source_line = {
                "id": midpoint_name + "__" + f"source_{idx}_{box_id}",
                "from_id": midpoint_name,
                "to_id": f"source_{idx}_{box_id}",
                "box_id": box_id,
                "type": 'transmission',
                "geometry": LineString([component_midpoints[idx]['geometry'], component_midpoints[idx]['geometry']]),  # line connecting to effective powerplant (line length = 0)
            }  # add details for effective powerplant
            box_simple_network_edges = box_simple_network_edges.append(midpoint_source_line, ignore_index=True)  # add effective powerplant connector
            box_simple_network_nodes = box_simple_network_nodes.append({'id':f'midpoint_source_{idx}_{box_id}', 'type':'source', 'box_id':box_id, 'capacity_mw':midpoint['total_mw'], 'geometry':midpoint['geometry']}, ignore_index=True)  # add effective powerplant


            # effective target
            midpoint_target_line = {
                "id": midpoint_name + "__" + f"target_{idx}_{box_id}",
                "from_id": midpoint_name,
                "to_id": f"target_{idx}_{box_id}",
                "box_id": box_id,
                "type": 'transmission',
                "geometry": LineString([component_midpoints[idx]['geometry'], component_midpoints[idx]['geometry']]),  # line connecting to effective target (line length = 0)
            }  # add details for effective target
            box_simple_network_edges = box_simple_network_edges.append(midpoint_target_line, ignore_index=True)  # add effective target connector
            box_simple_network_nodes = box_simple_network_nodes.append({'id':f'midpoint_target_{idx}_{box_id}', 'type':'target', 'box_id':box_id, 'gdp':midpoint['total_gdp'], 'geometry':midpoint['geometry']}, ignore_index=True)  # add effective target

        if len(box_simple_network_edges) == 0:  # if empty, write dummy
            cols_edges = ['id', 'box_id','geometry']
            box_simple_network_edges = gpd.GeoDataFrame(columns=cols_edges)
            box_simple_network_edges.loc[0, :] = [None] * len(cols_edges)
            box_simple_network_edges["box_id"] = box_id
        if len(box_simple_network_nodes) == 0:  # if empty, write dummy
            cols_nodes = cols_edges.copy()
            box_simple_network_nodes = gpd.GeoDataFrame(columns=cols_nodes)
            box_simple_network_nodes.loc[0, :] = [None] * len(cols_nodes)
            box_simple_network_nodes["box_id"] = box_id

        box_simple_network_edges.to_file(simple_file, layer='edges', driver='GPKG')
        box_simple_network_nodes.to_file(simple_file, layer='nodes', driver='GPKG')





