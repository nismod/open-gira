import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_make_exposure_tif():
    common.run_test(
        "make_exposure_tif",
        (
            "snakemake results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/raster/ "
            "-j1"
        ),
    )
