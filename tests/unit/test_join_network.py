import os
import sys

from . import runner


def test_join_network():
    runner.run_test(
        "join_network",
        (
            "snakemake results/djibouti-latest_filter-road/nodes.geoparquet "
            "results/djibouti-latest_filter-road/edges.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
