from . import runner


def test_join_data():
    runner.run_snakemake_test(
        "join_data",
        ("results/djibouti-latest_filter-road_hazard-aqueduct-river.geoparquet",),
    )
