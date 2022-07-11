import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_create_bbox_extracts():
    common.run_test(
        "create_bbox_extracts",
        "snakemake results/json/tanzania-mini_extracts.geojson -j1 --keep-target-files",
    )
