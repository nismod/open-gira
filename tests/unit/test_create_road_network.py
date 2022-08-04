import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_create_road_network():
    common.run_test(
        "create_road_network",
        (
            "snakemake results/geoparquet/djibouti-latest_filter-road/processed/slice-0_nodes.geoparquet "
            "results/geoparquet/djibouti-latest_filter-road/processed/slice-0_edges.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
