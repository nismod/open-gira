from . import runner


def test_estimate_wind_fields():
    runner.run_snakemake_test(
        "estimate_wind_fields",
        (
            "results/power/slice/1030/storms/IBTrACS/max_wind_field.nc",
        )
    )
