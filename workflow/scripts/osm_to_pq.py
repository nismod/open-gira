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

import geopandas
import pandas
import osmium
import shapely.geometry as shape
import shapely.ops as shape_ops


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
                    {"node": n, "point": shape.Point((n.lon, n.lat))}
                )
                if n.lon in node_index.keys():
                    node_index[n.lon][n.lat] = n
                else:
                    node_index[n.lon] = {n.lat: n}

        # Cut into bounding box multilinestring
        linestring = shape.linestring.LineString(locations)
        # constrain to bounding box
        lines = bbox.intersection(linestring)  # MULTILINESTRING | LINESTRING
        # split by shared nodes
        # GEOMETRYCOLLECTION of LINESTRINGs
        shared_points = shape.MultiPoint([n["point"] for n in shared_nodes_used])
        if lines.intersects(shared_points):
            segments = shape_ops.split(lines, shared_points)
        else:
            try:
                segments = shape.GeometryCollection(lines.geoms)
            except AttributeError:  # Single line
                segments = shape.GeometryCollection([lines])
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
                    "way_id": w.id,
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


if __name__ == "__main__":
    try:
        pbf_path = snakemake.input[0]
        output_path = snakemake.output[0]
        keep_tags = snakemake.config["keep_tags"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        pbf_path, output_path, keep_tags = sys.argv[1:]
        # pbf_path = '../../results/slices/tanzania-mini_filter-highway-core/slice-2.osm.pbf'
        # output_path = '../../results/test.geoparquet'
        # keep_tags = 'highway'

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    logging.info(f"Converting {pbf_path} to .geoparquet.")

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    box = osmium.io.Reader(pbf_path).header().box()
    bbox = shape.box(
        box.bottom_left.lon, box.bottom_left.lat, box.top_right.lon, box.top_right.lat
    )

    p = WayParser()
    shared_nodes = p.find_shared_nodes(pbf_path)

    h = WaySlicer(
        shared_nodes=shared_nodes,
        tags_to_preserve=keep_tags.replace(" ", "").split(","),
    )
    h.apply_file(pbf_path, locations=True)
    logging.info(
        f"Complete: {len(h.output_data)} segments from {len(Counter(w['way_id'] for w in h.output_data))} ways."
    )

    geopandas.GeoDataFrame(h.output_data).to_parquet(output_path)
