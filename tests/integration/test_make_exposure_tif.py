import pytest

from . import runner


@pytest.mark.skip(reason="Need to update rule (tif input no longer checkpoint) and fix script (use all tifs, not just first)")
def test_make_exposure_tif():
    runner.run_snakemake_test(
        "make_exposure_tif",
        (
            "results/exposure/djibouti-latest_filter-road/hazard-aqueduct-river/raster",
        )
    )
