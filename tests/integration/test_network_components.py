from . import runner


def test_network_components():
    runner.run_snakemake_test(
        "network_components",
        (
            "results/djibouti-latest_filter-road/component_population.svg",
            "results/djibouti-latest_filter-road/network_map_by_component.png",
            "results/djibouti-latest_filter-road/components.parquet",
        )
    )
