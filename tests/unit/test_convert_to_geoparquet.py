import os
import sys

from . import runner


def test_convert_to_geoparquet():
    runner.run_test(
        "convert_to_geoparquet",
        (
            "snakemake results/geoparquet/djibouti-latest_filter-road/raw/slice-0_edges.geoparquet "
            "results/geoparquet/djibouti-latest_filter-road/raw/slice-0_nodes.geoparquet "
            "-j1 --keep-target-files"
        ),
    )
