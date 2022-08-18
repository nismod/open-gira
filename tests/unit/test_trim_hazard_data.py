import os
import sys

from . import runner


def test_trim_hazard_data():
    runner.run_snakemake_test(
        "trim_hazard_data",
        (
            "results/input/hazard-aqueduct-river/djibouti-latest",
        )
    )
