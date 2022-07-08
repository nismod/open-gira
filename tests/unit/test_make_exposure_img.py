"""
These tests include downloads of the entire ne_50m admin boundary and ne_10m ocean maps.
It would be better if those maps were sliced down to just the bounding box covered by tanzania-mini.osm.pbf.

The image that is made does not contain coastline or admin boundary data, making it a poor test.
It would be better if the tests were overhauled to use a segment that _does_ have this information.
"""

import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_make_exposure_img():
    common.run_test(
        'make_exposure_img',
        (
            'snakemake results/exposure/tanzania-mini_filter-highway-core/hazard-aqueduct-river/img/ '
            '-j1'
        )
    )