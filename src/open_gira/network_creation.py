import logging
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from scipy.spatial import cKDTree
from shapely.geometry import LineString

import snkit


def create_network(
    edges: gpd.GeoDataFrame,
    nodes: gpd.GeoDataFrame = None,
    id_prefix: str = ""
) -> snkit.network.Network:
    """
    Create snkit network from edges and (optional) nodes and clean the result.

    Arguments:
        edges (gpd.GeoDataFrame): Expected to contain geometry column of linestrings
        nodes (gpd.GeoDataFrame): Optional nodes to include. snkit will try to snap to edges

    Returns:
        snkit.network.Network: Built network
    """

    logging.info("Starting network creation")

    # drop edges with no geometry
    empty_idx = edges.geometry.apply(lambda e: e is None or e.is_empty)
    if empty_idx.sum():
        empty_edges = edges[empty_idx]
        logging.info(f"Found {len(empty_edges)} empty edges.")
        logging.info(empty_edges)
        edges = edges[~empty_idx].copy()

    logging.info("Creating network")
    network = snkit.Network(nodes, edges)

    logging.info("Splitting multilines")
    network = snkit.network.split_multilinestrings(network)

    if nodes is not None and not nodes.empty:
        # check we have only point nodes
        assert set(network.nodes.geometry.type.values) == {"Point"}

        logging.info("Dropping duplicate geometries")
        # silence shapely.ops.split ShapelyDeprecationWarning regarding:
        # shapley.ops.split failure on split by empty geometry collection
        # this is currently caught by snkit.network.split_edge_at_points,
        # but won't be for shapely==2.0
        warnings.filterwarnings(
            "ignore",
            message=(
                ".*GeometryTypeError will derive from ShapelyError "
                "and not TypeError or ValueError in Shapely 2.0*"
            )
        )
        network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)

        logging.info("Snapping nodes to edges")
        network = snkit.network.snap_nodes(network)

    logging.info("Adding endpoints")
    network = snkit.network.add_endpoints(network)

    logging.info("Splitting edges at nodes")
    network = snkit.network.split_edges_at_nodes(network)

    # check we have only linestrings
    assert set(network.edges.geometry.type.values) == {"LineString"}

    logging.info("Renaming nodes and edges")
    network = snkit.network.add_ids(network, edge_prefix=id_prefix, node_prefix=id_prefix)

    logging.info("Creating network topology")
    network = snkit.network.add_topology(network, id_col="id")

    return network


def duplicate_reverse_and_append_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """
    Given edges with `from_id`, `to_id`, `from_iso_a4` and `to_iso_a3` columns,
    create duplicate edges with direction reversed.

    Args:
        edges: Table of edges to reverse and append to

    Returns:
        Table consisting of original edges and their reversed duplicates.
    """
    reversed_edges = edges.copy()
    reversed_edges.from_id = edges.to_id
    reversed_edges.to_id = edges.from_id
    reversed_edges.from_iso_a3 = edges.to_iso_a3
    reversed_edges.to_iso_a3 = edges.from_iso_a3
    return pd.concat([edges, reversed_edges])


def clean_maxspeed(value: str, default_km_h: float, min_km_h = 20, max_km_h = 140) -> float:
    """
    Cast, check and return the value of OSM maxspeed tag.

    Args:
        value: Speed limit value to clean
        default_km_h: Where data is missing or obviously wrong, return this value.

    Returns:
        Hopefully a sensible numeric speed limit value in km h-1
    """

    try:
        speed_km_h = float(value)
    except ValueError:
        return default_km_h

    if np.isnan(speed_km_h):
        return default_km_h
    elif (speed_km_h < min_km_h) or (speed_km_h > max_km_h):
        return default_km_h
    elif (speed_km_h >= min_km_h) and (speed_km_h <= max_km_h):
        return speed_km_h
    else:
        raise RuntimeError(f"Unforeseen consequences with {speed_km_h=}")


def preprocess_road_network(
    nodes_path: str,
    edges_path: str,
    filter_iso_a3: set[str],
    cost_USD_t_km: float,
    cost_USD_t_h: float,
    directional: float,
    default_max_speed_km_h: float,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Preprocess road network data into a suitable format for multi-modal routing.
    Find length of edges, calculate cost. Relabel IDs to include mode so they're unique across networks.

    Args:
        nodes_path: Path to nodes geoparquet file on disk
        edges_path: Path to edges geoparquet file on disk
        filter_iso_a3: ISO A3 codes of countries to keep features for
        cost_USD_t_km: Cost of transporting goods in USD per tonne km (broadly, cost of fuel)
        cost_USD_t_h: Cost of transporting goods in USD per tonne h (broadly, wages)
        directional: Whether to duplicate, reverse and append edges (from_id and to_id switched)
        default_max_speed_km_h: Value to gap fill speed limit data with

    Returns:
        Nodes and edges as two GeoDataFrames.
    """

    edges = gpd.read_parquet(edges_path)
    # at least one foot in the countries in question
    edges = edges[edges.from_iso_a3.isin(filter_iso_a3) | edges.to_iso_a3.isin(filter_iso_a3)]
    edges["distance_km"] = edges.geometry.to_crs(edges.estimate_utm_crs()).length / 1_000

    edges["mode"] = "road"
    edges["max_speed_km_h"] = edges.tag_maxspeed.apply(clean_maxspeed, args=(default_max_speed_km_h,))
    edges["avg_speed_km_h"] = edges.max_speed_km_h.apply(lambda x: np.clip(2/3 * x, None, default_max_speed_km_h))

    edges["cost_USD_t"] = cost_USD_t_km * edges["distance_km"] + cost_USD_t_h * edges["distance_km"] * 1 / edges["avg_speed_km_h"]
    edges["id"] = edges.apply(lambda row: f"{row['mode']}_{row['id']}", axis=1)
    edges["to_id"] = edges.apply(lambda row: f"{row['mode']}_{row['to_id']}", axis=1)
    edges["from_id"] = edges.apply(lambda row: f"{row['mode']}_{row['from_id']}", axis=1)

    if directional:
        edges = duplicate_reverse_and_append_edges(edges)

    nodes = gpd.read_parquet(nodes_path)
    nodes["mode"] = "road"
    nodes["id"] = nodes.apply(lambda row: f"{row['mode']}_{row['id']}", axis=1)

    return nodes, edges


def preprocess_rail_network(
    nodes_path: str,
    edges_path: str,
    filter_iso_a3: set[str],
    cost_USD_t_km: float,
    cost_USD_t_h: float,
    directional: float,
    avg_speed_km_h: float,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Preprocess rail network data into a suitable format for multi-modal routing.
    Find length of edges, calculate cost. Relabel IDs to include mode so they're unique across networks.

    Args:
        nodes_path: Path to nodes geoparquet file on disk
        edges_path: Path to edges geoparquet file on disk
        filter_iso_a3: ISO A3 codes of countries to keep features for
        cost_USD_t_km: Cost of transporting goods in USD per tonne km (broadly, cost of fuel)
        cost_USD_t_h: Cost of transporting goods in USD per tonne h (broadly, wages)
        directional: Whether to duplicate, reverse and append edges (from_id and to_id switched)
        avg_speed_km_h: Average speed of trains across network, used for all edges

    Returns:
        Nodes and edges as two GeoDataFrames.
    """

    edges = gpd.read_parquet(edges_path)
    # at least one foot in the countries in question
    edges = edges[edges.from_iso_a3.isin(filter_iso_a3) | edges.to_iso_a3.isin(filter_iso_a3)]
    edges["distance_km"] = edges.geometry.to_crs(edges.estimate_utm_crs()).length / 1_000

    edges["mode"] = "rail"

    edges["cost_USD_t"] = cost_USD_t_km * edges["distance_km"] + cost_USD_t_h * edges["distance_km"] * 1 / avg_speed_km_h
    edges["id"] = edges.apply(lambda row: f"{row['mode']}_{row['id']}", axis=1)
    edges["to_id"] = edges.apply(lambda row: f"{row['mode']}_{row['to_id']}", axis=1)
    edges["from_id"] = edges.apply(lambda row: f"{row['mode']}_{row['from_id']}", axis=1)

    if directional:
        edges = duplicate_reverse_and_append_edges(edges)

    nodes = gpd.read_parquet(nodes_path)
    nodes["mode"] = "rail"
    nodes["id"] = nodes.apply(lambda row: f"{row['mode']}_{row['id']}", axis=1)

    return nodes, edges


def preprocess_maritime_network(nodes_path: str, edges_path: str) -> tuple[gpd.GeoDataFrame, pd.DataFrame]:
    """
    Preprocess maritime network data into a suitable format for multi-modal routing.
    Relabel IDs so they're unique across networks.

    Args:
        nodes_path: Path to nodes geoparquet file on disk
        edges_path: Path to edges parquet file on disk

    Returns:
        Nodes GeoDataFrame and edges DataFrame.
    """
    edges = pd.read_parquet(edges_path).rename(columns={"from_iso3": "from_iso_a3", "to_iso3": "to_iso_a3"})
    edges["mode"] = "maritime"
    edges["cost_USD_t"] = edges["distance_km"] * edges["cost_USD_t_km"]

    nodes = gpd.read_parquet(nodes_path)
    nodes = nodes.rename(columns={"iso3": "iso_a3"})
    nodes = nodes.drop(columns=["Continent_Code"])
    ports_mask = nodes.infra == "port"

    # we want to connect our road and rail nodes to the port_land node of the port_in, port_out, port_land trifecta
    nodes.loc[ports_mask, "id"] = nodes.loc[ports_mask, :].apply(lambda row: f"{row.id}_land", axis=1)

    return nodes, edges


def find_nearest_points(
    a: gpd.GeoDataFrame,
    b: gpd.GeoDataFrame,
    b_id_col: str,
) -> gpd.GeoDataFrame:
    """
    Given two GeoDataFrames of point locations, `a` and `b`, for each point in `a`, find the closest in `b`.

    Modified from:
    https://gis.stackexchange.com/questions/222315/finding-nearest-point-in-other-geodataframe-using-geopandas

    a: Table of points to start from, we run over every row here
    b: Table of candidate closest points
    b_id_col: Name of column in b identifying points

    Returns:
        `a`, joined with the values from `b_id_col` from the points of `b` which are closest
    """

    # find nearest point in b for each and every point in a
    tree = cKDTree(b.geometry.get_coordinates().to_numpy())
    distances, indicies = tree.query(a.geometry.get_coordinates().to_numpy(), k=1)
    nearest_points = b.iloc[indicies][[b_id_col, "geometry"]] \
        .reset_index(drop=True).rename(columns={"geometry": "nearest_node_geometry"})

    return pd.concat(
        [
            a.reset_index(drop=True),
            nearest_points,
            pd.Series(data=distances, name='distance')
        ],
        axis=1
    )


def create_edges_to_nearest_nodes(
    a: gpd.GeoDataFrame,
    b: gpd.GeoDataFrame,
    max_distance_m: float,
    projected_coordinate_system: pyproj.crs.crs.CRS,
) -> gpd.GeoDataFrame:
    """
    Given two sets of nodes, a and b, loop through nodes in a, finding the closest
    node in b (which is less than `max_distance_m` away). Create a linear linestring
    connecting these points.

    Args:
        a: Table of nodes to connect from, containing GeoSeries of point locations.
            Must contain "id", "iso_a3" and "geometry" columns.
        b: Table of candidate notes to connect to, containing Geoseries of point locations.
            Must contain "id" and "geometry" columns.
        max_distance_m: Edges only created if their span in metres is equal to or less than this value.
        projected_coordinate_system: Project points to this CRS (must use metres!) before estimating distances.

    Returns:
        Table of linking edges.
    """

    point_pairs = find_nearest_points(
        a.to_crs(projected_coordinate_system),
        b.to_crs(projected_coordinate_system).rename(columns={"id": "nearest_node_id"}),
        "nearest_node_id"
    ).rename(columns={"distance": "distance_m"})

    point_pairs = point_pairs[point_pairs.distance_m < max_distance_m]

    edges = point_pairs.apply(
        lambda row: {
            "from_id": row.id,
            "to_id": row.nearest_node_id,
            "from_iso_a3": row.iso_a3,
            "to_iso_a3": row.iso_a3,  # assume link does not cross a border
            "geometry": LineString([row.geometry, row.nearest_node_geometry]),
            "distance_m": row.distance_m
        },
        axis=1,
        result_type="expand"
    )

    return gpd.GeoDataFrame(edges).reset_index(drop=True).set_crs(projected_coordinate_system)


def find_importing_node_id(row: pd.Series, exporting_country: str) -> str:
    """
    Return the node id lying in the importing country

    Args:
        row: Table row with columns from_iso_a3, to_iso_a3, from_id and to_id
        exporting_country: ISO A3 code of exporting country

    Returns
        node id
    """
    if row.from_iso_a3 == exporting_country and row.to_iso_a3 != exporting_country:
        return row.to_id
    elif row.from_iso_a3 != exporting_country and row.to_iso_a3 == exporting_country:
        return row.from_id
    else:
        raise RuntimeError


def create_edges_to_destination_countries(
    origin_nodes: gpd.GeoDataFrame,
    destination_country_nodes: gpd.GeoDataFrame,
    cost_USD_t: float = 1E6,
) -> gpd.GeoDataFrame:
    """
    Create edges between nodes within the same country of zero cost.

    Args:
        origin_nodes: Table of origin nodes to create edges from
        destination_country_nodes: Table of destination nodes to connect to, one per country
        cost_USD_t: Cost of traversing this edge, in USD per tonne. Defaults to a very high
            value, so that route allocations will only use these edges as a final connection
            to the destination, rather than a general means of traversal.

    Returns:
        Table of edges of same length as origin_nodes, connecting these to destination_country_nodes
    """

    assert len(destination_country_nodes.iso_a3) == len(destination_country_nodes.iso_a3.unique())
    assert origin_nodes.crs == destination_country_nodes.crs

    def make_edge(row: pd.Series) -> dict:
        destination = destination_country_nodes.set_index("iso_a3").loc[row.iso_a3]
        return {
            "from_id": row.id,
            "to_id": destination.id,
            "from_iso_a3": row.iso_a3,
            "to_iso_a3": row.iso_a3,
            "mode": "imaginary",
            "geometry": LineString(
                [
                    row.geometry,
                    destination.geometry
                ]
            ),
            "cost_USD_t": cost_USD_t
        }

    edges = origin_nodes.apply(
        make_edge,
        axis=1,
        result_type="expand"
    )

    return gpd.GeoDataFrame(edges).reset_index(drop=True).set_crs(origin_nodes.crs)


def path_edges_from_ordered_id_list(path_node_ids: list[str], edges: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    If `path_node_ids` are sequential nodes forming a path through a graph made of `edges`,
    return the subset of ordered edges connecting these nodes.

    Args:
        path_node_ids: Sequential node ids of some path through graph
        edges: Set of edges containing possible path edges

    Returns:
        Ordered subset of `edges` corresponding to path prescribed by `path_node_ids`
    """
    route_edges = []
    for i, from_node_id in enumerate(path_node_ids[:-1]):
        to_node_id = path_node_ids[i + 1]
        edge = edges[(edges.from_id == from_node_id) & (edges.to_id == to_node_id)]
        assert len(edge) == 1
        route_edges.append(edge)
    return pd.concat(route_edges)
