import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_join_data():
    # N.B. tests/config/config.yaml has an entry for slice_count
    # override it here to match the number of slices we have as input
    # this is necessary because the join_data rule is parameterised by
    # slice_count
    common.run_test(
        "join_data",
        (
            "snakemake results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet "
            "-j1 --keep-target-files --config slice_count=4"
        ),
    )
