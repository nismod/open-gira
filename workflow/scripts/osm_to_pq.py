"""
Convert OSM.pbf to parquet using osmium, geopandas

Take an osm.pbf file and decompose it into a list of road segments.
Each segment is defined by a start node and an end node,
and consists of a section of road that has no junctions except at the
nodes.

The process is as follows:
* Make a list of all of the node references in an OSM slice
* TODO: Create nodes for the intersections of way and OSM slice bounding box and mark as shared
* Drop unique node references
* TODO: Drop nodes that lie outside OSM bounding box
* For each way, break it into segments that join its start, shared, or end nodes
* Convert each segment to a LineString and save as .geoparquet

"""
import logging
import math
from collections import Counter

import geopandas
import osmium
import shapely.geometry.linestring


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

    def get_key_nodes(self):
        return self.shared_nodes.keys()

    def way(self, w):
        # Find segments
        shared_node_set = {w.nodes[0].ref, *self.shared_nodes.keys(), w.nodes[-1].ref}
        # print(f"way={w.id}; shared node count={len([n for n in w.nodes if n.ref in shared_node_set])}")
        segment = []
        segments = {}
        segment_id = 0
        for node in w.nodes:
            if node.ref in shared_node_set and len(segment) >= 2:
                # end segment
                segment.append(node)
                segments[segment_id] = segment
                segment_id += 1
                segment = [node]
            else:
                # start/middle
                segment.append(node)

        # Prepare information for all segments
        base_input = {}
        for k in self.tags_to_preserve:
            base_input[f"tag_{k}"] = w.tags[k] if k in w.tags else None

        # Calculate segment output
        for s_id, segment in segments.items():
            """
            osmium will not let us create a linestring from a mutable object (not sure why!).
            That means we have to create the linestring manually.
            """
            locations = [(n.lon, n.lat) for n in segment]
            line = shapely.geometry.linestring.LineString(locations)
            self.output_data.append({
                'geometry': line,
                'way_id': w.id,
                'segment_id': s_id,
                **base_input,
                **self.shared_node_to_dict(segment[0], 'start_node_'),
                **self.shared_node_to_dict(segment[-1], 'end_node_')
            })

    def shared_node_to_dict(self, shared_node, prefix):
        """
        Return a dictionary of node information with entries prefixed by prefix
        Parameters
        ----------
        shared_node :Node: node to process
        prefix :str: prefix to use for dictionary entries

        Returns
        -------
        dictionary of node reference, longitude, latitude, and degree
        """
        try:
            degree = self.shared_nodes[shared_node.ref]
        except KeyError:
            degree = 1
        return {
            f'{prefix}reference': shared_node.ref,
            f'{prefix}longitude': shared_node.lon,
            f'{prefix}latitude': shared_node.lat,
            f'{prefix}degree': degree
        }


if __name__ == '__main__':
    try:
        pbf_path = snakemake.input[0]
        output_path = snakemake.output[0]
        keep_tags = snakemake.config['keep_tags']
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        pbf_path, output_path, keep_tags = sys.argv[1:]

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    logging.info(f"Converting {pbf_path} to .geoparquet.")

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

    p = WayParser()
    shared_nodes = p.find_shared_nodes(pbf_path)

    h = WaySlicer(
        shared_nodes=shared_nodes,
        tags_to_preserve=keep_tags.replace(' ', '').split(',')
    )
    h.apply_file(pbf_path, locations=True)
    logging.info(
        f"Complete: {len(h.output_data)} segments from {len(Counter(w['way_id'] for w in h.output_data))} ways."
    )

    geopandas.GeoDataFrame(h.output_data).to_parquet(output_path)
