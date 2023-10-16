from . import runner


def test_estimate_wind_fields():
    runner.run_snakemake_test(
        "estimate_wind_fields",
        (
            "results/power/by_country/PRI/storms/IBTrACS/0/max_wind_field.nc",
        )
    )
