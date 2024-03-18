from . import runner


def test_filter_osm_data():
    runner.run_snakemake_test(
        "filter_osm_data",
        (
            "results/input/OSM/djibouti-latest_filter-road-secondary.osm.pbf",
        )
    )
