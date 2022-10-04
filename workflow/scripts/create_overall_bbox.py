import sys
import subprocess
import json
import re
import os.path

try:
    osm_file = snakemake.input[0]  # type: ignore
    results_dir = snakemake.config["output_dir"]  # type: ignore
    out_file = snakemake.output[0]  # type: ignore
except NameError:
    if len(sys.argv) != 4:
        raise RuntimeError(
            "Incorrect number of input args, 3 required. Args: .osm.pbf file, results directory, new .json file"
        )
    osm_file = sys.argv[1]
    results_dir = sys.argv[2]
    out_file = sys.argv[3]

bboxes = subprocess.check_output(["osmium", "fileinfo", osm_file, "-g", "header.boxes"])

# TODO: support multiple bounding boxes (osmium help says these print a multiline output to the above command)
box = re.match("^\\((-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+)", bboxes.decode())
if box:
    content = {
        "extracts": [
            {
                "bbox": [float(x) for x in list(box.groups())],
                "output": os.path.basename(osm_file),
            }
        ],
    }
    with open(out_file, "w") as out:
        out.write(json.dumps(content))
