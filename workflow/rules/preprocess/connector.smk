"""
Finds connector points for all boxes
"""
from open_gira.process_power_functions import adj

def adjacent(wildcards):
    output_dir = wildcards.OUTPUT_DIR
    with open(f"{output_dir}/processed/world_boxes_metadata.json") as filejson:
        world_boxes_metadata = json.load(filejson)

    num_cols = world_boxes_metadata["num_cols"]
    tot_boxes = world_boxes_metadata["tot_boxes"]

    # Adjacent box numbers
    box_id = int(wildcards.BOX)
    neighbour_box_ids = adj(box_id, num_cols, tot_boxes)
    edges = [
        f"{output_dir}/processed/power/{b}/edges_{b}.parquet"
        for b in neighbour_box_ids
    ]
    nodes = [
        f"{output_dir}/processed/power/{b}/nodes_{b}.parquet"
        for b in neighbour_box_ids
    ]
    return edges + nodes

rule process_connector:
    conda: "../../../environment.yml"
    input:
        adjacent,
        edges="{OUTPUT_DIR}/processed/power/{BOX}/edges_{BOX}.parquet",
        nodes="{OUTPUT_DIR}/processed/power/{BOX}/nodes_{BOX}.parquet",
        global_metadata=rules.world_splitter.output.global_metadata,
    output:
        connector="{OUTPUT_DIR}/processed/power/{BOX}/connector_{BOX}.json",
    script:
        "../../scripts/process/process_power_5_connector.py"
