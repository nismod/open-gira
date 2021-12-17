import os
from collections import defaultdict
from glob import glob

import geopandas as gpd
import networkx as nx
import pandas as pd
from tqdm.notebook import tqdm
from importing_modules import *
from process_power_functions import *

fname = os.path.join("data","processed","world_network.gpkg")


#changedir()  # todo remove

def read_network(fname):
    # read edges
    edges = gpd.read_file(fname, layer='edges')

    # create unique link id from from/to node ids
    edges['link'] = edges.apply(lambda e: "__".join(sorted([e.from_id, e.to_id])), axis=1)

    # read nodes
    nodes = gpd.read_file(fname, layer='nodes')

    # add gdp to target nodes
    targets = gpd.read_file(fname.replace('network','targets'))
    nodes = nodes.merge(targets[['id', 'gdp']], on='id', how='left')

    # add capacity to plant/generation/source nodes
    plants = gpd.read_file(fname.replace('network','plants'))
    if nodes['capacity_mw'].all() == plants['capacity_mw'].all():
        nodes = nodes.drop(columns="capacity_mw")  # prevent merge issue
    nodes = nodes.merge(plants[['id', 'capacity_mw']], on='id', how='left')

    return nodes, edges


def create_graph(nodes, edges):
    G = nx.Graph()
    G.add_nodes_from(
        (
            n.id,
            {"id": n.id, "type": n.type, "gdp": n.gdp, "capacity_mw": n.capacity_mw}
        )
        for n in nodes.itertuples()
    )
    G.add_edges_from(
        (
            e.from_id,
            e.to_id,
            {"id": e.id, "length_m": e.geometry.length}
        )
        for e in edges.itertuples()
    )
    return G


def assign_node_edge_gdp(G):
    components = list(nx.connected_components(G))
    node_gdps = []
    edge_gdps = []

    for component in tqdm(components):
        c_node_gdp, c_edge_gdp = assign_component_gdp(G, component)
        if len(c_node_gdp):
            node_gdps.append(c_node_gdp)
        if len(c_edge_gdp):
            edge_gdps.append(c_edge_gdp)

    node_gdp = pd.concat(node_gdps)
    edge_gdp = pd.concat(edge_gdps)

    return node_gdp, edge_gdp


def edge_link_ids_from_nodes(G, route_nodes):
    next_nodes = iter(route_nodes)
    next(next_nodes)
    return [
        "__".join(sorted([u, v]))
        for u, v in zip(route_nodes, next_nodes)
    ]


def assign_component_gdp(G, component):
    # create connected component subgraph
    c = G.subgraph(component).copy()
    # nodes dataframe
    c_nodes = pd.DataFrame(n for _, n in c.nodes(data=True))

    # component total gdp and capacity
    c_gdp = c_nodes.gdp.sum()
    c_cap = c_nodes.capacity_mw.sum()

    # assign GDP to nodes in proportion to capacity
    def assign_source_gdp(n):
        if n.type == 'source':
            return n.capacity_mw * c_gdp / c_cap
        return n.gdp

    c_nodes['gdp'] = c_nodes.apply(assign_source_gdp, axis=1)

    # assign GDP "flow" along shortest path from source to target, sharing source
    # gdp proportionally between targets
    edge_links = defaultdict(int)
    sources = c_nodes[c_nodes.type == 'source'].copy()
    targets = c_nodes[c_nodes.type == 'target'].copy()
    if len(sources) and len(targets):
        for u in tqdm(sources.itertuples(), total=len(sources)):
            paths = nx.shortest_path(c, source=u.id)
            for v in targets.itertuples():
                path = paths[v.id]
                path_gdp = u.gdp * (v.gdp / c_gdp)
                for link_id in edge_link_ids_from_nodes(c, path):
                    edge_links[link_id] += path_gdp

    c_edges = pd.DataFrame({'link': k, 'gdp': v} for k, v in edge_links.items())

    return c_nodes[['id', 'gdp']], c_edges


#%% run

print("importing")
nodes, edges = read_network(fname)
timer(start)

print("creating graph")
G = create_graph(nodes, edges)
timer(start)

print("assigning node edges gdp")
node_gdp, edge_gdp = assign_node_edge_gdp(G)
timer(start)

out_fname = fname.replace('network', 'network_with_gdp')

print("node processing")
nodes = nodes.drop('gdp', axis=1).merge(node_gdp, on='id')
timer(start)

print("writing nodes to file")
nodes.to_file(
    out_fname,
    layer='nodes',
    driver='GPKG')
timer(start)

print("merging edges")
edges = edges.merge(edge_gdp, on='link')
timer(start)

print("writing edges to file")
edges.to_file(
    fname.replace('network', 'network_with_gdp'),
    layer='edges',
    driver='GPKG')

end = time.time()
print(f"\nTime to run file: {(end - start)/60:0.2f} minutes")