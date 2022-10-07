from . import runner


def test_network_raster():
    runner.run_snakemake_test(
        "network_raster",
        (
            "results/splits/djibouti-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet",
        )
    )
