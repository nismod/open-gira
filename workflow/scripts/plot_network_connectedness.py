"""
Plotting tools for analysing network connectedness
"""

import sys
from typing import Iterable, Tuple

import pandas as pd
import geopandas as gpd
import networkx
import numpy as np
import matplotlib.pyplot as plt


def network_components(edges: pd.DataFrame, nodes: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, list[set]]:
    """
    Label nodes and edges with a component_id, where for any given component_id, all nodes
    labelled such are connected to one another.
    """

    # build graph
    graph = networkx.Graph()
    graph.add_nodes_from(n.node_id for n in nodes.itertuples())
    graph.add_edges_from((e.from_node_id, e.to_node_id) for e in edges.itertuples())

    # list of sets of nodes
    # each set is a component, an 'island' of nodes
    components = list(networkx.connected_components(graph))

    # assign component_id to edges and nodes
    component_id_col: str = "component_id"
    for index, component in enumerate(components):
        # indexing the edges with both the from and the to nodes should be redundant here
        edges.loc[(edges.from_node_id.isin(component) | edges.to_node_id.isin(component)), component_id_col] = index
        nodes.loc[nodes["node_id"].isin(component), component_id_col] = index

    return edges, nodes, components


def plot_component_size(components: Iterable[set[str]]) -> Tuple[plt.Figure, plt.axis]:
    """Plot a bar chart of the size of each network component."""

    f, ax = plt.subplots(figsize=(8, 6))
    component_sizes = sorted([len(n) for n in components])
    ax.bar(range(len(component_sizes)), component_sizes)
    ax.set_yscale('log')
    ax.set_ylabel('Nodes in component')
    ax.set_title('Size of components in network')

    # annotate plot with some stats on how many nodes in x of y components
    thresholds = np.round(np.logspace(0, 1, 4))
    for i, n_large_components in enumerate(thresholds[thresholds < len(components)].astype(int)):

        largest_components = component_sizes[-n_large_components:]
        pct_nodes_largest_comp = 100 * sum(largest_components) / sum(component_sizes)
        ax.text(
            0.06,
            0.9 - 0.05 * i,
            f"Top {n_large_components} of {len(components)} components contain "
            f"{sum(largest_components)} of {sum(component_sizes)} nodes ({pct_nodes_largest_comp:.2f}%)",
            transform=ax.transAxes
        )

    return f, ax


if __name__ == "__main__":

    try:
        nodes_path = snakemake.input["nodes"]
        edges_path = snakemake.input["edges"]
        plot_path = snakemake.output["plot"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        nodes_path, edges_path, plot_path = sys.argv[1:]
        # nodes_path = ../../results/tanzania-mini_filter-highway-core/road_edges.geoparquet
        # edges_path = ../../results/tanzania-mini_filter-highway-core/road_edges.geoparquet
        # plot_path = ../../results/tanzania-mini_filter-highway-core/road_connectedness.pdf

    edges = gpd.read_parquet(edges_path)
    nodes = gpd.read_parquet(nodes_path)
    edges, nodes, components = network_components(edges, nodes)
    f, ax = plot_component_size(components)
    f.savefig(plot_path)
