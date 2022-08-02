import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_create_overall_bbox():
    common.run_test(
        "create_overall_bbox",
        "snakemake results/json/djibouti-latest.json -j1 --keep-target-files",
    )
