"""
Plotting tools for analysing network connectedness
"""

import os
import re
import sys
from collections import defaultdict
from typing import Iterable

import datashader as ds
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import snkit
import spatialpandas
from datashader.utils import export_image
from matplotlib import colors

from open_gira.utils import natural_sort


def random_cmap(n) -> colors.LinearSegmentedColormap:
    """Create a colormap with n random colors."""

    rgb = np.random.rand(n, 3)
    # bolt on a column of ones for the alpha channel
    rgba = np.c_[rgb, np.ones(n)]
    return colors.LinearSegmentedColormap.from_list("random", rgba)


def plot_components_map(network: snkit.network.Network) -> ds.transfer_functions.Image:
    """Draw a map of all the edges in the network. Uses datashader for speed."""

    # use spatialpandas rather than geopandas ;) lol thanks python ecosystem
    # datashader requires input data to be of this type
    edges = spatialpandas.GeoDataFrame(network.edges)

    # this will probably break near the poles or with the whole planet,
    # but try and estimate the canvas size from the geometry bounding box
    pixels_per_degree = 100
    x0, y0, x1, y1 = edges.geometry.total_bounds
    cvs = ds.Canvas(
        plot_width=int(pixels_per_degree * (x1 - x0)),
        plot_height=int(pixels_per_degree * (y1 - y0)),
    )

    # datashader is all about aggregating data down to a raster...
    # take the mean value of component_id for all the edges on that pixel
    # unfortunately the mode is not (yet) available, but this works fine for now
    agg = cvs.line(edges, geometry="geometry", agg=ds.mean("component_id"))

    # color the map according to a random colormap (should show each island in the network)
    cmap = random_cmap(len(edges))
    # if we ony have one component, colour the edges white
    # N.B. datashader will fail to shade with a conventional colourmap in this situation
    if len(set(edges.component_id.values)) == 1:
        cmap = ["white"]

    return ds.transfer_functions.shade(agg, how="log", cmap=cmap)


def plot_component_size(components: Iterable[set[str]]) -> tuple[plt.Figure, plt.axis]:
    """Plot a bar chart of the size of each network component."""

    f, ax = plt.subplots(figsize=(8, 6))
    component_sizes = sorted([len(n) for n in components])
    ax.bar(range(len(component_sizes)), component_sizes)
    ax.set_yscale("log")
    ax.set_ylabel("Nodes in component")
    ax.set_title("Size of components in network")

    # annotate plot with some stats on how many nodes in x of y components
    thresholds = np.round(np.logspace(0, 1, 4))
    for i, n_large_components in enumerate(
        thresholds[thresholds < len(components)].astype(int)
    ):

        largest_components = component_sizes[-n_large_components:]
        pct_nodes_largest_comp = 100 * sum(largest_components) / sum(component_sizes)
        ax.text(
            0.06,
            0.9 - 0.05 * i,
            f"Top {n_large_components} of {len(components)} components contain "
            f"{sum(largest_components)} of {sum(component_sizes)} nodes ({pct_nodes_largest_comp:.2f}%)",
            transform=ax.transAxes,
        )

    return f, ax


if __name__ == "__main__":

    try:
        nodes_path = snakemake.input["nodes"]  # type: ignore
        edges_path = snakemake.input["edges"]  # type: ignore
        population_plot_path = snakemake.output["component_population"]  # type: ignore
        map_path = snakemake.output["component_map"]  # type: ignore
        components_path = snakemake.output["component_data"]  # type: ignore
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        (
            nodes_path,
            edges_path,
            population_plot_path,
            map_path,
            components_path,
        ) = sys.argv[1:]
        # nodes_path = ../../results/tanzania-mini_filter-road/road_edges.geoparquet
        # edges_path = ../../results/tanzania-mini_filter-road/road_edges.geoparquet
        # population_plot_path = ../../results/tanzania-mini_filter-road/road_component_population.pdf
        # map_path = ../../results/tanzania-mini_filter-road/road_network_map_by_component.pdf
        # components_path = ../../results/tanzania-mini_filter-road/road_components.parquet

    # build a network from files on disk
    network = snkit.network.Network(
        edges=gpd.read_parquet(edges_path),
        nodes=gpd.read_parquet(nodes_path)
    )

    # extract the component data
    # 'edge_ids' or 'node_ids' -> component_id -> set of element ids
    component_map: dict[int, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

    # N.B. sort the entries to make testing easier
    for component_id in sorted(set(network.edges.component_id)):
        component_mask = network.edges.component_id == component_id
        component_map["edge_ids"][component_id] = natural_sort(network.edges.id[component_mask])
        node_ids = set(network.edges.from_id[component_mask]) | set(network.edges.to_id[component_mask])
        component_map["node_ids"][component_id] = natural_sort(node_ids)

    # build a pandas dataframe from the components
    df = pd.DataFrame(component_map)
    df.index = df.index.rename("component_id")

    # write out nodes and edges of each component
    df.to_parquet(components_path)

    # bar chart of component edge populations
    edges_by_component: list[set[str]] = [c["edges"] for c in component_map.values()]
    fig, ax = plot_component_size(df.edge_ids)
    fig.savefig(population_plot_path)

    # map coloured by component id
    map_image = plot_components_map(network)
    stem, ext = os.path.splitext(map_path)
    export_image(map_image, stem, fmt=ext, background="black")
