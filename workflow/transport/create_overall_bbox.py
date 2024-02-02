import json
import os.path
import re
import subprocess


# do not permit values outside this range
# otherwise, we create many empty slices for the polar regions
MIN_LAT = -60
MAX_LAT = 72

osm_file = snakemake.input.osm_pbf
results_dir = snakemake.wildcards.OUTPUT_DIR
out_file = snakemake.output.bbox

bboxes = subprocess.check_output(["osmium", "fileinfo", osm_file, "-g", "header.boxes"])

# TODO: support multiple bounding boxes (osmium help says these print a multiline output to the above command)
box = re.match("^\\((-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+)", bboxes.decode())

if box:
    min_long, min_lat, max_long, max_lat = [float(x) for x in list(box.groups())]
    content = {
        "extracts": [
            {
                "bbox": [min_long, max(min_lat, MIN_LAT), max_long, min(max_lat, MAX_LAT)],
                "output": os.path.basename(osm_file),
            }
        ],
    }
    with open(out_file, "w") as out:
        out.write(json.dumps(content))
