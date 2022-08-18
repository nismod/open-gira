"""
These tests include downloads of the entire ne_50m admin boundary and ne_10m ocean maps.
It would be better if those maps were sliced down to just the bounding box covered by djibouti-latest.osm.pbf.
"""

import os
import sys

from . import runner


def test_make_exposure_img():
    runner.run_snakemake_test(
        "make_exposure_img",
        (
            "results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/img",
        )
    )
