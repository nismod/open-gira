from . import runner


def test_download_hazard_datasets():
    runner.run_snakemake_test(
        "download_hazard_datasets", ("results/input/hazard-aqueduct-river/raw",)
    )
