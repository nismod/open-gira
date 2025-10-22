from . import runner


def test_rasterise_osm_network():
    runner.run_snakemake_test(
        "rasterise_osm_network",
        (
            "results/splits/djibouti-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet",
        ),
    )
