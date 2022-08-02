import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_join_data():
    common.run_test(
        "join_data",
        (
            "snakemake results/djibouti-latest_filter-road_hazard-aqueduct-river.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
