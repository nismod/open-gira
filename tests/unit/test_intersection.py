import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_intersection():
    common.run_test(
        "intersection",
        (
            "snakemake results/splits/djibouti-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
