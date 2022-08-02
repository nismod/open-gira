import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_join_network():
    common.run_test(
        "join_network",
        (
            "snakemake results/djibouti-latest_filter-road/nodes.geoparquet "
            "results/djibouti-latest_filter-road/edges.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
