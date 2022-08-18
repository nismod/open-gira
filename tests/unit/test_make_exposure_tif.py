from . import runner


def test_make_exposure_tif():
    runner.run_snakemake_test(
        "make_exposure_tif",
        (
            "results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/raster",
        )
    )
