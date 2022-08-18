import os
import sys

from . import runner


def test_create_bbox_extracts():
    runner.run_snakemake_test(
        "create_bbox_extracts",
        (
            "results/json/djibouti-latest_extracts.geojson",
        )
    )
