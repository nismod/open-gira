import os
import sys

from . import runner


def test_intersection():
    runner.run_snakemake_test(
        "intersection",
        (
            "results/splits/djibouti-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet",
        )
    )
