from . import runner


def test_electricity_grid_damages():
    runner.run_snakemake_test(
        "electricity_grid_damages",
        (
            "results/power/by_country/PRI/exposure/IBTrACS/0/",
            "results/power/by_country/PRI/disruption/IBTrACS/0/",
        )
    )
