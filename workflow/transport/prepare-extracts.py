"""
Take bounding box defined in input GeoJSON file and split into slices, write
slices to list of output files.
"""

import json
import math
import os
import sys
import re
from os.path import splitext as splitext_orig

from open_gira.utils import natural_sort


def splitext(string, exts=None):
    dataset, ext = splitext_orig(string)
    if ext == "":
        return dataset, exts
    exts = ext if exts is None else ext + exts
    return splitext(dataset, exts)


def slice_subextracts(initial_bbox, n):
    """
    Arguments
    ---------
    initial_bbox: List of xmin, xmax, ymin, ymax (list[float])
    n: int

    Returns:
    --------
    List of bounding boxes, one for each slice (list[list[float]])
    """

    xmin0, ymin0, xmax0, ymax0 = tuple(initial_bbox)
    dx = (xmax0 - xmin0) / n
    dy = (ymax0 - ymin0) / n

    bboxes = []
    for ix in range(n):
        xmin = xmin0 + ix * dx
        xmax = xmin0 + (ix + 1) * dx
        for iy in range(n):
            ymin = ymin0 + iy * dy
            ymax = ymin0 + (iy + 1) * dy

            bboxes.append([xmin, ymin, xmax, ymax])

    return bboxes


try:
    slice_count = int(snakemake.config["slice_count"])  # type: ignore
    input_json_path = snakemake.input[0]  # type: ignore
    extract_paths: list[str] = snakemake.output  # type: ignore
except NameError:
    slice_count = int(sys.argv[1])
    input_json_path = sys.argv[2]
    extract_paths: list[str] = sys.argv[3:]

# sort for reproducibility
extract_paths = natural_sort(extract_paths)

# check inputs
# fail if we have more than one out dir in the supplied paths
out_dir, = set(map(os.path.dirname, extract_paths))
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# check we've got a square number
# N.B. this is also checked in the Snakefile
n = math.sqrt(slice_count)
if n % 1 > 0 and n != 1:
    raise ValueError("Total slice count must be a square number or 1.")
elif slice_count == 1:
    n = slice_count
else:
    n = int(n)

# read the initial bbox
with open(input_json_path, "r") as fp:
    # fail if more than one initial bounding box
    extract, = json.load(fp)["extracts"]

# generate all slice bboxes given sqrt(slice_count)
slice_bboxes: list[list[float]] = slice_subextracts(extract["bbox"], n)

# write out slice sub-extracts for the paths we've been given
for extract_path in extract_paths:
    extract_filename = os.path.basename(extract_path)
    i, = re.match(r"slice-(\d+)", extract_filename).groups()

    slice_json = {
        "directory": ".",
        "extracts": [
            {
                "bbox": slice_bboxes[int(i)],
                "output": extract_filename.replace(".geojson", ".osm.pbf")
            }
        ]
    }

    with open(extract_path, "w") as fp:
        json.dump(slice_json, fp, indent=4)
