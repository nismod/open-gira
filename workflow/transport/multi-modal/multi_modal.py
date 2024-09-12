
import geopandas as gpd
import numpy as np
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from open_gira.network_creation import (
    duplicate_reverse_and_append_edges, preprocess_road_network,
    preprocess_rail_network, create_edges_to_nearest_nodes,
    find_importing_node_id, create_edges_to_destination_countries
)
from open_gira.routing import DESTINATION_LINK_COST_USD_T

matplotlib.use("Agg")

if __name__ == "__main__":

    study_country: str = snakemake.config["study_country_iso_a3"]

    print("Preprocessing road network...")
    road_nodes, road_edges = preprocess_road_network(
        snakemake.input.road_network_nodes,
        snakemake.input.road_network_edges,
        {study_country,},
        snakemake.config["road_cost_USD_t_km"],
        snakemake.config["road_cost_USD_t_h"],
        True,
        snakemake.config["road_default_speed_limit_km_h"]
    )

    print("Preprocessing rail network...")
    rail_nodes, rail_edges = preprocess_rail_network(
        snakemake.input.rail_network_nodes,
        snakemake.input.rail_network_edges,
        {study_country,},
        snakemake.config["rail_cost_USD_t_km"],
        snakemake.config["rail_cost_USD_t_h"],
        True,
        snakemake.config["rail_average_freight_speed_km_h"]
    )

    print("Reading maritime network...")
    maritime_nodes = gpd.read_parquet(snakemake.input.maritime_nodes) 
    maritime_edges = gpd.read_parquet(snakemake.input.maritime_edges)

    maximum_intermodal_connection_metres = 2_000

    print("Making intermodal connections...")
    # road-rail
    rail_road_edges = create_edges_to_nearest_nodes(
        rail_nodes.loc[rail_nodes.station == True, ["id", "iso_a3", "geometry"]],
        road_nodes.loc[:, ["id", "geometry"]],
        maximum_intermodal_connection_metres,
        rail_nodes.estimate_utm_crs()
    ).to_crs(epsg=4326)
    rail_road_edges["mode"] = "road_rail"

    # road-maritime
    maritime_road_edges = create_edges_to_nearest_nodes(
        maritime_nodes.loc[
            (maritime_nodes.infra == "port") & (maritime_nodes.iso_a3 == study_country),
            ["id", "iso_a3", "geometry"]
        ],
        road_nodes.loc[:, ["id", "geometry"]],
        maximum_intermodal_connection_metres,
        road_nodes.estimate_utm_crs()
    ).to_crs(epsg=4326)
    maritime_road_edges["mode"] = "maritime_road"

    # rail-maritime
    maritime_rail_edges = create_edges_to_nearest_nodes(
        maritime_nodes.loc[
            (maritime_nodes.infra == "port") & (maritime_nodes.iso_a3 == study_country),
            ["id", "iso_a3", "geometry"]
        ],
        rail_nodes.loc[rail_nodes.station == True, ["id", "geometry"]],
        maximum_intermodal_connection_metres,
        road_nodes.estimate_utm_crs()
    ).to_crs(epsg=4326)
    maritime_rail_edges["mode"] = "maritime_rail"

    intermodal_edges = pd.concat(
        [
            rail_road_edges,
            maritime_road_edges,
            maritime_rail_edges
        ]
    ).to_crs(epsg=4326)
    # as the maritime edges are directional, we're making road, rail and intermodal directional too (so duplicate)
    intermodal_edges = duplicate_reverse_and_append_edges(intermodal_edges)

    intermodal_cost_USD_t = snakemake.config["intermodal_cost_USD_t"]
    intermodal_edges["cost_USD_t"] = intermodal_edges["mode"].map(intermodal_cost_USD_t)

    # concatenate different kinds of nodes and edges
    node_cols = ["id", "iso_a3", "geometry"]
    nodes = gpd.GeoDataFrame(
        pd.concat(
            [
                road_nodes.loc[:, node_cols],
                rail_nodes.loc[:, node_cols],
                maritime_nodes.loc[:, node_cols]
            ]
        ),
        crs=4326
    )

    edge_cols = ["from_id", "to_id", "from_iso_a3", "to_iso_a3", "mode", "cost_USD_t", "geometry"]
    edges = pd.concat(
        [
            intermodal_edges.loc[:, edge_cols],
            road_edges.loc[:, edge_cols],
            rail_edges.loc[:, edge_cols],
            maritime_edges.loc[:, edge_cols]
        ]
    )
    # add nodes for destination countries (not null, not origin country)
    # neighbouring countries will have destination node connected to border crossings
    countries = nodes.iso_a3.unique()
    countries = countries[countries != np.array(None)]
    countries = set(countries)
    countries.remove(study_country)

    admin_boundaries = gpd.read_parquet(snakemake.input.admin_boundaries)
    country_nodes = admin_boundaries.set_index("GID_0").loc[list(countries), ["geometry"]] \
        .sort_index().reset_index().rename(columns={"GID_0": "iso_a3"})
    country_nodes["id"] = country_nodes.apply(lambda row: f"GID_0_{row.iso_a3}", axis=1)
    country_nodes.geometry = country_nodes.geometry.representative_point()
    destination_country_nodes = country_nodes.loc[:, ["id", "iso_a3", "geometry"]]

    nodes = pd.concat(
        [
            nodes,
            destination_country_nodes
        ]
    )

    # find nodes which lie on far side of border crossing
    border_crossing_mask = \
        (edges.from_iso_a3 != edges.to_iso_a3) \
        & ((edges.from_iso_a3 == study_country) | (edges.to_iso_a3 == study_country)) \
        & ((edges["mode"] == "road") | (edges["mode"] == "rail")) \
        
    importing_node_ids = edges[border_crossing_mask].apply(find_importing_node_id, exporting_country=study_country, axis=1)
    importing_nodes = nodes.set_index("id").loc[importing_node_ids].reset_index()
    # two importing nodes are labelled as THA, drop these
    importing_nodes = importing_nodes[importing_nodes.iso_a3 != study_country]

    # connect these nodes to their containing country
    land_border_to_importing_country_edges = \
        create_edges_to_destination_countries(importing_nodes, destination_country_nodes)
    land_border_to_importing_country_edges.plot()

    print("Plotting land border crossings for inspection...")
    # plot THA land border crossing points for sanity
    to_plot = importing_nodes
    country_ints, labels = pd.factorize(to_plot["iso_a3"])
    unique_country_ints = []
    cmap = plt.get_cmap("viridis")
    colours = [cmap(x) for x in country_ints / max(country_ints)]
    [unique_country_ints.append(c) for c in colours if c not in unique_country_ints]
    colour_map = dict(zip(labels, unique_country_ints))
    f, ax = plt.subplots()
    to_plot.plot(color=colours, ax=ax)
    ax.set_title("Land border crossing points")
    patches = [mpatches.Patch(color=colour, label=label) for label, colour in colour_map.items()]
    ax.legend(handles=patches)
    f.savefig(snakemake.output.border_crossing_plot)

    print("Making terminal connections to destination countries...")
    # connect foreign ports to their country with new edges
    foreign_ports = maritime_nodes[maritime_nodes.infra=="port"]
    foreign_ports = foreign_ports[foreign_ports.iso_a3 != study_country]
    port_to_importing_countries_edges = create_edges_to_destination_countries(
        foreign_ports,
        destination_country_nodes,
        DESTINATION_LINK_COST_USD_T
    )

    # add in edges connecting destination countries to THA land borders and foreign ports
    edges = pd.concat(
        [
            edges.loc[:, edge_cols],
            duplicate_reverse_and_append_edges(land_border_to_importing_country_edges.loc[:, edge_cols]),
            duplicate_reverse_and_append_edges(port_to_importing_countries_edges.loc[:, edge_cols]),
        ]
    ).reset_index(drop=True)

    # there are duplicate edges (repeated from_id -> to_id pairs), drop these here
    edges["unique_edge_id"] = edges.apply(lambda row: f"{row.from_id}_{row.to_id}", axis=1)
    edges = edges[~edges.unique_edge_id.duplicated(keep="first")].drop(columns=["unique_edge_id"])

    print("Write out network to disk as geoparquet...")
    # write out global multi-modal transport network to disk
    # reset indicies to 0-start integers
    # these will correspond to igraph's internal edge/vertex ids
    nodes = nodes.reset_index(drop=True)
    nodes.to_parquet(snakemake.output.nodes)
    edges = edges.reset_index(drop=True)
    edges.to_parquet(snakemake.output.edges)