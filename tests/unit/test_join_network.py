from . import runner


def test_join_network():
    runner.run_snakemake_test(
        "join_network",
        (
            "results/djibouti-latest_filter-road/nodes.geoparquet",
            "results/djibouti-latest_filter-road/edges.geoparquet",
        )
    )
