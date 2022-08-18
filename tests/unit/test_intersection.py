import os
import sys

from . import runner


def test_intersection():
    runner.run_test(
        "intersection",
        (
            "snakemake results/splits/djibouti-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
