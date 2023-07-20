"""
Convert OSM.pbf to parquet using osmium, geopandas

Take an osm.pbf file and decompose it into a list of road segments.
Each segment is defined by a start node and an end node,
and consists of a section of road that has no junctions except at the
nodes.

The process is as follows:
* Make a list of all node references in an OSM slice
* Convert to shapely LINESTRING
* Intersect with bounding box
* Split by any nodes shared with other ways
* Save with way and node details
"""

import logging
import sys
from collections import Counter
from typing import Union

import geopandas
import osmium
import pandas
import shapely.ops as shape_ops
from shapely.errors import GeometryTypeError
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
from shapely.geometry import box as shapely_box
from shapely.geometry.collection import GeometryCollection


class WayParser(osmium.SimpleHandler):
    """
    Generate a list of node references for identifying non-unique nodes
    """

    def __init__(self):
        osmium.SimpleHandler.__init__(self)
        self.node_list = []
        self.shared_nodes = None

    def way(self, w):
        [self.node_list.append(n.ref) for n in w.nodes]

    def find_shared_nodes(self, file_path, locations=False, **kwargs):
        if self.shared_nodes is not None:
            return self.shared_nodes
        self.shared_nodes = {}
        self.apply_file(file_path, locations=locations, **kwargs)
        node_counts = Counter(self.node_list)
        for k, v in node_counts.items():
            if v > 1:
                self.shared_nodes[k] = v
        return self.shared_nodes


class NodeParser(osmium.SimpleHandler):
    """
    Extract nodes from OSM data
    """

    def __init__(self, tags_to_preserve):
        osmium.SimpleHandler.__init__(self)
        self.output_data = []
        self.tags_to_preserve = tags_to_preserve

    def node(self, n):
        """
        Process an individual node and add it to the output list
        """

        base_input = {}
        for k in self.tags_to_preserve:
            base_input[f"tag_{k}"] = n.tags[k] if k in n.tags else None

        # create shapely geometry from osm object
        # https://docs.osmcode.org/pyosmium/latest/ref_osm.html#osmium.osm.Location
        point = Point(n.location.lon, n.location.lat)

        self.output_data.append(
            {
                "geometry": point,
                "osm_node_id": n.id,
                **base_input,
            }
        )


class WaySlicer(osmium.SimpleHandler):
    """
    Slice up ways into segments by shared nodes
    @param List<int> shared_nodes - list of nodes that are shared with other ways in the network
    @param List<string> tags_to_preserve - list of osmium tags to keep in the output
    """

    def __init__(self, shared_nodes, tags_to_preserve):
        osmium.SimpleHandler.__init__(self)
        self.output_data = []
        self.shared_nodes = shared_nodes
        self.tags_to_preserve = tags_to_preserve

    def way(self, w):
        if len(w.nodes) < 2:
            # not enough points in this way to create a linestring, short circuit
            return

        # Prepare information for all segments
        base_input = {}
        for k in self.tags_to_preserve:
            base_input[f"tag_{k}"] = w.tags[k] if k in w.tags else None

        # Parse nodes for relevant properties (one pass)
        locations = []
        shared_nodes_used = []
        node_index = {}
        for n in w.nodes:
            locations.append((n.lon, n.lat))
            if n.ref in self.shared_nodes.keys():
                shared_nodes_used.append(
                    {"node": n, "point": Point((n.lon, n.lat))}
                )
                if n.lon in node_index.keys():
                    node_index[n.lon][n.lat] = n
                else:
                    node_index[n.lon] = {n.lat: n}

        # generate linestring and constrain to bounding box
        cut_to_bbox: Union[Point, Linestring, MultiLineString] = bbox.intersection(LineString(locations))
        if not isinstance(cut_to_bbox, (LineString, MultiLineString)):
            # the intersection can return a point when the way has one end
            # outside the bbox, and its other end on the bbox edge
            # in this case, discard the way
            return

        # split by shared nodes
        # GEOMETRYCOLLECTION of LINESTRINGS
        shared_points = MultiPoint([n["point"] for n in shared_nodes_used])
        if cut_to_bbox.intersects(shared_points):
            segments = shape_ops.split(cut_to_bbox, shared_points)
        else:
            try:
                segments = GeometryCollection(cut_to_bbox.geoms)
            except AttributeError:  # Single line
                segments = GeometryCollection([cut_to_bbox])
        s_id = 0
        for line in segments.geoms:
            # Determine start and end nodes from shared nodes or invent if bbox clipping
            prefixes = ["start_node_", "end_node_"]
            nodes = []
            for i in range(2):
                n = line.coords[i]
                try:
                    node = self.get_node_by_coords(n, prefixes[i], node_index)
                except KeyError:
                    # Not found, must be a node we created by clipping to bbox
                    # TODO: investigate possible bug: these nodes not just at bbox!
                    node = {
                        f"{prefixes[i]}reference": pandas.NA,
                        f"{prefixes[i]}longitude": n[0],
                        f"{prefixes[i]}latitude": n[1],
                        f"{prefixes[i]}degree": 1,
                    }
                nodes.append(node)

            self.output_data.append(
                {
                    "geometry": line,
                    "osm_way_id": w.id,
                    "segment_id": s_id,
                    **base_input,
                    **nodes[0],
                    **nodes[1],
                }
            )
            s_id += 1

    def get_node_by_coords(self, coords, prefix, node_list):
        """
        Return a dictionary of node information with entries prefixed by prefix
        Parameters
        ----------
        coords :Node: node to process
        prefix :str: prefix to use for dictionary entries
        node_list :Dict<Node>: dict of candidate nodes that shared_node might match, indexed by lon, lat

        Returns
        -------
        dictionary of node reference, longitude, latitude, and degree
        """
        node = node_list[coords[0]][coords[1]]  # KeyError is caught by parent function
        if (
            node.ref not in self.shared_nodes.keys()
        ):  # The shared_nodes should all have degree > 1
            raise RuntimeError(f"Node {node.ref} not found in shared_nodes keys.")
        degree = self.shared_nodes[node.ref]
        return {
            f"{prefix}reference": node.ref,
            f"{prefix}longitude": node.lon,
            f"{prefix}latitude": node.lat,
            f"{prefix}degree": degree,
        }


def empty_gdf() -> geopandas.GeoDataFrame:
    """
    Create an return an empty GeoDataFrame. Must explicitly specify columns
    (despite empty list) to permit saving as geoparquet.
    """
    return geopandas.GeoDataFrame([], columns=["geometry"])


if __name__ == "__main__":
    try:
        pbf_path = snakemake.input["pbf"]  # type: ignore
        edges_path = snakemake.output["edges"]  # type: ignore
        nodes_path = snakemake.output["nodes"]  # type: ignore
        keep_tags = snakemake.params["keep_tags"]  # type: ignore
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        pbf_path, edges_path, nodes_path, keep_tags = sys.argv[1:]
        # pbf_path = 'results/slices/tanzania-mini_filter-road/slice-2.osm.pbf'
        # edges_path = 'results/slice-2.geoparquet'
        # nodes_path = 'results/slice-2.geoparquet'
        # keep_tags = 'highway, railway'

        # process comma separated string into list of strings
        keep_tags: list = keep_tags.replace(" ", "").split(",")

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    logging.info(f"Converting {pbf_path} to .geoparquet.")

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    box = osmium.io.Reader(pbf_path).header().box()
    bbox = shapely_box(
        box.bottom_left.lon, box.bottom_left.lat, box.top_right.lon, box.top_right.lat
    )

    p = WayParser()
    shared_nodes = p.find_shared_nodes(pbf_path)

    h = WaySlicer(
        shared_nodes=shared_nodes,
        tags_to_preserve=keep_tags,
    )
    h.apply_file(pbf_path, locations=True)

    if len(h.output_data) != 0:
        edges = geopandas.GeoDataFrame(h.output_data)
        edges = edges.set_crs(epsg=4326)
    else:
        edges = empty_gdf()
    logging.info(
        f"Complete: {len(h.output_data)} segments from {len(Counter(w['osm_way_id'] for w in h.output_data))} ways."
    )

    n = NodeParser(
        tags_to_preserve=keep_tags,
    )
    n.apply_file(pbf_path, locations=True)
    if len(n.output_data) != 0:
        nodes = geopandas.GeoDataFrame(n.output_data)
        nodes = nodes.set_crs(epsg=4326)
    else:
        nodes = empty_gdf()

    logging.info(f"Complete: {len(n.output_data)} nodes.")

    # write to disk -- even if empty
    edges.to_parquet(edges_path)
    nodes.to_parquet(nodes_path)
