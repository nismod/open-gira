from . import runner


def test_event_set_direct_damages():
    runner.run_snakemake_test(
        "event_set_direct_damages",
        (
            "results/direct_damages/djibouti-latest_filter-road/hazard-jba-event/fraction/slice-0.gpq",
            "results/direct_damages/djibouti-latest_filter-road/hazard-jba-event/cost/slice-0.gpq"
        )
    )
