"""Extracts the box_ids within the selected region

"""
import os


## Process Config Inputs ##
REGIONS = config["regions"]
if len(REGIONS) == 0:
    print("Inputting all regions")
    REGIONS = ["EP", "NA", "NI", "SI", "SP", "WP"]


SAMPLES = list(range(config["sample_upper"] + 1))
if config["samples_indiv"] != "None":
    print("Using specified samples")
    SAMPLES = config["samples_indiv"]
if len(SAMPLES) == 0:
    print("Samples incorrectly specified")

STORMS = config["specific_storm_analysis"]
if STORMS == 'None':
    STORMS = None



region_box = expand(
    os.path.join("data", "intersection", "regions", "{region}_boxes.txt"),
    region=REGIONS,
)


rule intersect_region_boxes:
    input:
        os.path.join(config['output_dir'], "power_processed", "world_boxes_metadata.txt"),
        #os.path.join(config['output_dir'], "power_processed", "world_boxes.gpkg"),  # removed because opening on QGIS tampers with metadata
        os.path.join(
            config['output_dir'], "input", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"
        ),
        expand(
            os.path.join(
                config['output_dir'],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "gridfinder_{box_id}.gpkg",
            ),
            box_id=all_boxes,
        ),
    output:
        os.path.join(config['output_dir'], "power_intersection", "regions", "{region}_boxes.txt"),
    params:
        region="{region}",
        output_dir = config['output_dir']
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_1_regions.py")


rule intersect_regions:
    input:
        region_box,
