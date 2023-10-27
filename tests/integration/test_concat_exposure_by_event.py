from . import runner


def test_concat_exposure_by_event():
    runner.run_snakemake_test(
        "concat_exposure_by_event",
        (
            "results/power/by_country/PRI/exposure/IBTrACS/exposure_by_event.parquet",
        )
    )
