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
from collections import Counter, defaultdict
from typing import Union

import geopandas
import osmium
import pandas
import shapely
import shapely.ops as shape_ops
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

    def way(self, way: osmium.osm.Way) -> None:
        """
        Add the node references from `way` to `self.node_list`. This method is
        called from `self.apply_file` for each way in the file.

        Args:
            way: Way to store node references for.
        """
        [self.node_list.append(node.ref) for node in way.nodes]

    def find_shared_nodes(
        self, file_path: str, locations: bool = False, **kwargs
    ) -> dict[int, int]:
        """
        Determine which nodes have a degree of 2 or greater.

        Args:
            file_path: Path to .osm.pbf file.
            locations: If locations is true, then a location handler will be
                used in `self.apply_file`, returning the node positions.

        Returns:
            Mapping from node ID to node degree for all 'shared_nodes'.
        """

        if self.shared_nodes is not None:
            return self.shared_nodes

        self.shared_nodes = {}
        self.apply_file(file_path, locations=locations, **kwargs)

        node_counts = Counter(self.node_list)
        for node_id, degree in node_counts.items():
            if degree > 1:
                self.shared_nodes[node_id] = degree

        return self.shared_nodes


class NodeParser(osmium.SimpleHandler):
    """
    Extract nodes from OSM data.
    """

    def __init__(self, tags_to_preserve: tuple[str]):
        osmium.SimpleHandler.__init__(self)
        self.output_data = []
        self.tags_to_preserve = tags_to_preserve

    def node(self, node: osmium.osm.Node) -> None:
        """
        Process an individual node and add it to the output list. This method
        is called from `self.apply_file` for each node in the file.

        Args:
            node: Node to process.
        """

        base_input = {}
        for tag_name in self.tags_to_preserve:
            base_input[f"tag_{tag_name}"] = (
                node.tags[tag_name] if tag_name in node.tags else None
            )

        # https://docs.osmcode.org/pyosmium/latest/ref_osm.html#osmium.osm.Location
        point: shapely.Point = Point(node.location.lon, node.location.lat)

        self.output_data.append(
            {
                "geometry": point,
                "osm_node_id": node.id,
                **base_input,
            }
        )


class WaySlicer(osmium.SimpleHandler):
    """
    Slice up OSM ways into segments, split at junctions. Clip ways to bounding box.

    Args:
        shared_nodes: Nodes with degree greater than or equal to two.
        tags_to_preserve: Way attributes we wish to retain.
        bounding_box: Bounding box to slice ways to. Any way extending beyond
            the bounding box will be lost from this slice.
    """

    def __init__(
        self,
        shared_nodes: dict[int, int],
        tags_to_preserve: tuple[str],
        bounding_box: shapely.Polygon,
    ):
        osmium.SimpleHandler.__init__(self)
        self.output_data = []
        self.shared_nodes = shared_nodes
        self.tags_to_preserve = tags_to_preserve
        self.bounding_box = bounding_box

    def way(self, way: osmium.osm.Way) -> None:  # noqa: C901
        """
        Slice given way into segments, splitting at junctions. Append processed
        way to `self.output_data`.

        Note that this method is called from `self.apply_file`.

        Args:
            way: Way to process.
        """

        if len(way.nodes) < 2:
            # not enough points in this way to create a linestring, short circuit
            return

        # extract tags for way
        base_input = {}
        for tag_name in self.tags_to_preserve:
            base_input[f"tag_{tag_name}"] = (
                way.tags[tag_name] if tag_name in way.tags else None
            )

        node_locations: tuple[float, float] = []
        shared_nodes_used: list[dict] = []
        # create parallel data structure to enable lookup by location
        nodes_indexed_by_location: defaultdict[
            float, dict[float, osmium.osm.NodeRef]
        ] = defaultdict(dict)

        # step through all nodes in way, check if any are shared with other ways
        for node in way.nodes:
            node_locations.append((node.lon, node.lat))

            # add it to the nodes dict indexed by longitude and latitude
            nodes_indexed_by_location[node.lon][node.lat] = node

            # node is shared between this way and another neighbouring way
            if node.ref in self.shared_nodes.keys():
                # add it to the 'used' shared nodes dict
                shared_nodes_used.append(
                    {"node": node, "point": Point((node.lon, node.lat))}
                )

        # Generate linestring and constrain to bounding box
        way_linestring_cut_to_bbox: Union[Point, LineString, MultiLineString] = (
            self.bounding_box.intersection(LineString(node_locations))
        )

        if way_linestring_cut_to_bbox.is_empty or (
            not isinstance(way_linestring_cut_to_bbox, (LineString, MultiLineString))
        ):
            # the intersection can return a point when the way has one end
            # outside the bbox, and its other end on the bbox edge
            # There are also degenerate ways with two co-located nodes which
            # are returned from intersection as an empty LineString
            # In either case, discard the way
            return

        # split at shared nodes, so we always have edges ending at junctions
        junction_points = MultiPoint([n["point"] for n in shared_nodes_used])

        if way_linestring_cut_to_bbox.intersects(junction_points):
            way_segments: shapely.GeometryCollection = shape_ops.split(
                way_linestring_cut_to_bbox, junction_points
            )
        else:
            try:
                way_segments = GeometryCollection(way_linestring_cut_to_bbox.geoms)
            except AttributeError:  # Single line
                way_segments = GeometryCollection([way_linestring_cut_to_bbox])

        # loop through segments in way
        # note that many ways will only have a single segment
        for segment_id, segment in enumerate(way_segments.geoms):
            # determine start and end nodes from shared nodes or invent if bbox has clipped way
            termini = []
            for terminus_index, prefix in ((0, "start_node_"), (-1, "end_node_")):
                longitude, latitude = segment.coords[terminus_index]

                try:
                    # an end of this segment is an existing shared node
                    terminus_data = self.get_node_by_coords(
                        longitude, latitude, prefix, nodes_indexed_by_location
                    )
                except KeyError:
                    # no record of a node at this location, must be a way clipped by the bounding box
                    logging.debug(
                        f"Clipped way {way.id} to bounding box at ({longitude:.2f}, {latitude:.2f})"
                    )
                    terminus_data = {
                        f"{prefix}reference": pandas.NA,
                        f"{prefix}longitude": longitude,
                        f"{prefix}latitude": latitude,
                        f"{prefix}degree": 1,
                    }

                termini.append(terminus_data)

            start, end = termini
            self.output_data.append(
                {
                    "geometry": segment,
                    "osm_way_id": way.id,
                    "segment_id": segment_id,
                    **base_input,
                    **start,
                    **end,
                }
            )

    def get_node_by_coords(
        self,
        longitude: float,
        latitude: float,
        prefix: str,
        nodes_indexed_by_location: dict[float, dict[float, osmium.osm.NodeRef]],
    ) -> dict:
        """
        Create a dictionary of node attributes, where keys are prefixed by `prefix`.

        Args:
            longitude: Node longitude.
            latitude: Node latitude.
            prefix: Prefix to use for dictionary keys.
            nodes_indexed_by_location: Node objects, indexed by longitude, latitude.

        Raises:
            KeyError: If `nodes_indexed_by_location` does not contain node at
                `longitude` and `latitude`.

        Returns:
            Node reference, longitude, latitude, and degree.
        """

        # note that we're indexing with floats here, which _should_ still be safe
        # given that they're immutable and we don't perform any operations on the floats
        node: osmium.osm.NodeRef = nodes_indexed_by_location[longitude][latitude]

        try:
            node_degree: int = self.shared_nodes[node.ref]
        except KeyError:
            # nodes that aren't shared should all have degree 1
            node_degree: int = 1

        return {
            f"{prefix}reference": node.ref,
            f"{prefix}longitude": node.lon,
            f"{prefix}latitude": node.lat,
            f"{prefix}degree": node_degree,
        }


def empty_gdf() -> geopandas.GeoDataFrame:
    """
    Create an return an empty GeoDataFrame. Must explicitly specify columns
    (despite empty list) to permit saving as geoparquet.
    """
    return geopandas.GeoDataFrame([], columns=["geometry"])


if __name__ == "__main__":
    pbf_path: str = snakemake.input["pbf"]  # noqa: F821
    edges_path: str = snakemake.output["edges"]  # noqa: F821
    nodes_path: str = snakemake.output["nodes"]  # noqa: F821
    keep_tags: tuple[str] = tuple(snakemake.params["keep_tags"])  # noqa: F821

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.DEBUG
    )
    logging.info(f"Converting {pbf_path} to .geoparquet.")

    # Ignore geopandas parquet implementation warnings
    import warnings

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    box = osmium.io.Reader(pbf_path).header().box()
    bbox = shapely_box(
        box.bottom_left.lon, box.bottom_left.lat, box.top_right.lon, box.top_right.lat
    )

    # EDGES

    # node ID -> node degree
    shared_nodes: dict[int, int] = WayParser().find_shared_nodes(pbf_path)

    way_slicer = WaySlicer(
        shared_nodes=shared_nodes, tags_to_preserve=keep_tags, bounding_box=bbox
    )
    way_slicer.apply_file(pbf_path, locations=True)

    if len(way_slicer.output_data) != 0:
        edges = geopandas.GeoDataFrame(way_slicer.output_data)
        edges = edges.set_crs(epsg=4326)
    else:
        edges = empty_gdf()

    logging.info(
        f"Complete: {len(way_slicer.output_data)} segments from "
        f"{len(Counter(way['osm_way_id'] for way in way_slicer.output_data))} ways."
    )

    # NODES

    node_parser = NodeParser(tags_to_preserve=keep_tags)
    node_parser.apply_file(pbf_path, locations=True)

    if len(node_parser.output_data) != 0:
        nodes = geopandas.GeoDataFrame(node_parser.output_data)
        nodes = nodes.set_crs(epsg=4326)

        # some nodes belong to ways which extended beyond the bounding box
        # we have sliced the ends off these ways, and now we discard their nodes
        nodes = nodes[nodes.intersects(bbox)]
    else:
        nodes = empty_gdf()

    logging.info(f"Complete: {len(node_parser.output_data)} nodes.")

    # write to disk, empty or not
    edges.to_parquet(edges_path)
    nodes.to_parquet(nodes_path)
