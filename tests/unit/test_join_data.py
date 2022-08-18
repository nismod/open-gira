import os
import sys

from . import runner


def test_join_data():
    runner.run_test(
        "join_data",
        (
            "snakemake results/djibouti-latest_filter-road_hazard-aqueduct-river.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
