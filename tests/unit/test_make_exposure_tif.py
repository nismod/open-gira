import os
import sys

from . import runner


def test_make_exposure_tif():
    runner.run_test(
        "make_exposure_tif",
        (
            "snakemake results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/raster/ "
            "-j1"
        ),
    )
