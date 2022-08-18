import os
import sys

from . import runner


def test_filter_osm_data():
    runner.run_test(
        "filter_osm_data",
        "snakemake results/input/djibouti-latest_filter-road.osm.pbf -j1 --keep-target-files",
    )
