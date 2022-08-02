import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_convert_to_geoparquet():
    common.run_test(
        "convert_to_geoparquet",
        (
            "snakemake results/geoparquet/djibouti-latest_filter-road/raw/slice-0_edges.geoparquet "
            "results/geoparquet/djibouti-latest_filter-road/raw/slice-0_nodes.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
