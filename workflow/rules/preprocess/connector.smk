"""
Finds connector points for all boxes
"""
from open_gira.process_power_functions import adj

def adjacent_box_nodes(wildcards):
    output_dir = wildcards.OUTPUT_DIR
    with open(f"{output_dir}/power/world_boxes_metadata.json") as filejson:
        world_boxes_metadata = json.load(filejson)

    num_cols = world_boxes_metadata["num_cols"]
    tot_boxes = world_boxes_metadata["tot_boxes"]

    # Adjacent box numbers
    box_id = int(wildcards.BOX)
    neighbour_box_ids = adj(box_id, num_cols, tot_boxes)
    nodes = [
        f"{output_dir}/power/slice/{b}/nodes_{b}.geoparquet"
        for b in neighbour_box_ids
    ]
    return nodes


rule process_connector:
    conda: "../../../environment.yml"
    input:
        adjacent_nodes=adjacent_box_nodes,
        edges="{OUTPUT_DIR}/power/slice/{BOX}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/slice/{BOX}/network/nodes.geoparquet",
        global_metadata=rules.world_splitter.output.global_metadata,
    output:
        connector="{OUTPUT_DIR}/power/slice/{BOX}/network/connector.json",
    script:
        "../../scripts/process/process_power_5_connector.py"
