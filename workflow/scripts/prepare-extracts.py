"""Slice collection of geometries definded in geoJSON file and write
slices in geoJSONfile.
"""

import json
import math
from os.path import join
from os.path import splitext as splitext_orig
import sys


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
    Iterator on bounding boxes, one for each slice (list[list[float]])
    """
    xmin0, ymin0, xmax0, ymax0 = tuple(initial_bbox)
    dx = (xmax0 - xmin0) / n
    dy = (ymax0 - ymin0) / n
    for ix in range(n):
        xmin = xmin0 + ix * dx
        xmax = xmin0 + (ix + 1) * dx
        for iy in range(n):
            ymin = ymin0 + iy * dy
            ymax = ymin0 + (iy + 1) * dy

            yield [xmin, ymin, xmax, ymax]


slice_count = int(sys.argv[2])
n = math.sqrt(slice_count)
if n % 1 > 0:
    raise ValueError('Total slice count must be a square number.')
else:
    n = int(n)

with open(sys.argv[1], "r") as fp:
    originaljsonfile = json.load(fp)

output_dir = originaljsonfile["directory"]
subextractsjson = {"directory": output_dir, "extracts": []}
for extract in originaljsonfile["extracts"]:
    dataset, ext = splitext(extract["output"])
    for n, bbox in enumerate(slice_subextracts(extract["bbox"], n)):
        subextractsjson["extracts"].append(
            {"bbox": bbox, "output": f"{dataset}-slice{n}{ext}"}
        )
    with open(join(sys.argv[3], dataset + "-extracts.geojson"), "w") as fp:
        json.dump(subextractsjson, fp, indent=4)
    subextractsjson["extracts"].clear()
