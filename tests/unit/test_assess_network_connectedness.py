from . import runner


def test_assess_network_connectedness():
    runner.run_snakemake_test(
        "assess_network_connectedness",
        (
            "results/djibouti-latest_filter-road/component_population.svg",
            "results/djibouti-latest_filter-road/network_map_by_component.png",
            "results/djibouti-latest_filter-road/components.parquet",
        )
    )
