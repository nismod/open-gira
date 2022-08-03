"""
These tests include downloads of the entire ne_50m admin boundary and ne_10m ocean maps.
It would be better if those maps were sliced down to just the bounding box covered by djibouti-latest.osm.pbf.
"""

import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_make_exposure_img():
    common.run_test(
        "make_exposure_img",
        (
            "snakemake results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/img/ "
            "-j1"
        ),
    )
