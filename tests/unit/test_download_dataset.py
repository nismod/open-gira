import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_download_dataset():
    common.run_test(
        "download_dataset",
        "snakemake results/input/tanzania-mini.osm.pbf -F -j1 --keep-target-files",
    )
