import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_filter_osm_data():
    common.run_test(
        "filter_osm_data",
        "snakemake results/input/tanzania-mini_filter-road.osm.pbf -j1 --keep-target-files",
    )
