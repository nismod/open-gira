from . import runner


def test_exposure_by_admin_region():
    runner.run_snakemake_test(
        "exposure_by_admin_region",
        (
            "results/power/by_country/PRI/exposure/IBTrACS/admin-level-1.geoparquet",
        )
    )
