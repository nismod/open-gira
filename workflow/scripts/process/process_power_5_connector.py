"""Documents the connections between boxes"""
import json
import os

import geopandas
from shapely.geometry import LineString

from open_gira.process_power_functions import adj


output_dir = snakemake.config["output_dir"]  # type: ignore
global_metadata_path = snakemake.input.global_metadata  # type: ignore
nodes_path = snakemake.input.nodes  # type: ignore
edges_path = snakemake.input.edges  # type: ignore
connector_path = snakemake.output.connector  # type: ignore
box_id = int(snakemake.wildcards.BOX)  # type: ignore


with open(global_metadata_path) as filejson:
    world_boxes_metadata = json.load(filejson)
num_cols = world_boxes_metadata["num_cols"]
tot_boxes = world_boxes_metadata["tot_boxes"]

# Adjacent box numbers
neighbour_box_ids = adj(box_id, num_cols, tot_boxes)

portal_dict = {}  # master dict
nodes = geopandas.read_parquet(nodes_path)

for neighbour_box_id in neighbour_box_ids:
    portal_lst_ex = []  # neighbour list
    neighbour_node_path = os.path.join(
        output_dir,
        "processed",
        "power",
        f"{neighbour_box_id}",
        f"network_box_{neighbour_box_id}.parquet",
    )
    if os.path.exists(neighbour_node_path):
        neighbour_nodes = geopandas.read_parquet(neighbour_node_path)
        duplicates = neighbour_nodes.merge(
            nodes,
            left_on="geometry",
            right_on="geometry",
            suffixes=("_ex", "_idx"),
        )
        update_lst = [
            {
                "edge_id": f"edge_X_{from_box}__{to_box}",
                "source_id": from_source_id,
                "link": f"cb_{from_box}__{to_box}",
                "type": "transmission",
                "from_id": from_idx,  # connection in box_idx
                "to_id": to_ex,  # connection in box next to box_idx
                "geometry": str(LineString([fromto_geom, fromto_geom])),
            }
            for from_idx, to_ex, from_box, to_box, from_source_id, fromto_geom in zip(
                duplicates.id_idx,
                duplicates.id_ex,
                duplicates.box_id_idx,
                duplicates.box_id_ex,
                duplicates.source_id_idx,
                duplicates.geometry,
            )
        ]

        # update master dict with neighbour dict
        portal_dict.update(
            {f"box_{neighbour_box_id}": update_lst}
        )

with open(connector_path, "w") as fp:
    json.dump(portal_dict, fp, indent=2)
