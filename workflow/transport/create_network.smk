"""
Generic network creation rules.
"""


import logging

import geopandas as gpd
import pandas as pd

import snkit


rule create_transport_network:
    """
    Take .geoparquet OSM files and output files of cleaned network nodes and edges
    """
    input:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_edges.geoparquet",
        admin="{OUTPUT_DIR}/input/admin-boundaries/gadm36_levels.gpkg",
    output:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_edges.geoparquet"
    params:
        # determine the network type from the filter, e.g. road, rail
        # example FILTER_SLUG values might be 'filter-road-tertiary' or 'filter-rail'
        network_type=lambda wildcards: wildcards.FILTER_SLUG.split('-')[1],
        # pass in the slice number so we can label edges and nodes with their slice
        # edge and node IDs should be unique across all slices
        slice_number=lambda wildcards: int(wildcards.SLICE_SLUG.replace('slice-', ''))
    script:
        # template the path string with a value from params (can't execute .replace in `script` context)
        "./create_{params.network_type}_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-road/processed/slice-0_edges.geoparquet
"""


def transport_network_paths_from_file(wildcards) -> list[str]:
    """
    Lookup composite network components from file and return list of paths to
    their edge files.
    """
    df = pd.read_csv(
        config["composite_network"][wildcards.COMPOSITE],
        comment="#",
    )
    edge_paths = df.apply(
        lambda row:
        f"{wildcards.OUTPUT_DIR}/{row.infrastructure_dataset}_filter-{row.network_filter}/edges.gpq",
        axis=1
    )
    node_paths = [path.replace("edges.gpq", "nodes.gpq") for path in edge_paths]
    return {"component_nodes": node_paths, "component_edges": edge_paths}


rule create_composite_transport_network:
    """
    Stitch together a set of transport networks into one file. These should be
    the same transport mode. May be used for creating networks of spatially
    varying density.
    """
    input:
        unpack(transport_network_paths_from_file)
    output:
        composite_nodes = "{OUTPUT_DIR}/composite_network/{COMPOSITE}/nodes.gpq",
        composite_edges = "{OUTPUT_DIR}/composite_network/{COMPOSITE}/edges.gpq"
    run:
        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Concatenate nodes and edges")
        nodes = []
        for node_path in input.component_nodes:
            nodes.append(gpd.read_parquet(node_path))
        edges = []
        for edge_path in input.component_edges:
            edges.append(gpd.read_parquet(edge_path))

        network = snkit.network.Network(
            nodes=pd.concat(nodes).reset_index(drop=True),
            edges=pd.concat(edges).reset_index(drop=True),
        )

        logging.info("Labelling edge ends with from/to node ids")
        network = snkit.network.add_topology(network)

        logging.info("Labelling edges and nodes with network component ids")
        network = snkit.network.add_component_ids(network)

        logging.info("Writing network to disk")
        network.nodes.to_parquet(output.composite_nodes)
        network.edges.to_parquet(output.composite_edges)

"""
Test with:
snakemake -c1 -- results/composite_network/south-east-asia-road/edges.gpq
"""
