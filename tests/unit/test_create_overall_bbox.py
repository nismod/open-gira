from . import runner


def test_create_overall_bbox():
    runner.run_snakemake_test(
        "create_overall_bbox",
        (
            "results/json/djibouti-latest.json",
        )
    )
