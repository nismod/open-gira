from . import runner


def test_electricity_grid_damages():
    runner.run_snakemake_test(
        "electricity_grid_damages",
        (
            "results/power/slice/1030/exposure/IBTrACS.nc",
        )
    )
