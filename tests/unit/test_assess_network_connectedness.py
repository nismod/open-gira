import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_assess_network_connectedness():
    common.run_test(
        "assess_network_connectedness",
        (
            "snakemake results/djibouti-latest_filter-road/component_population.svg "
            "results/djibouti-latest_filter-road/network_map_by_component.png "
            "results/djibouti-latest_filter-road/components.parquet "
            "-j1 --keep-target-files"
        ),
    )
