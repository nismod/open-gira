from . import runner


def test_aggregate_exposure_within_sample():
    runner.run_snakemake_test(
        "aggregate_exposure_within_sample",
        (
            "results/power/by_country/PRI/exposure/IBTrACS/0_length_m_by_event.pq",
            "results/power/by_country/PRI/exposure/IBTrACS/0_length_m_by_edge.pq",
        )
    )
