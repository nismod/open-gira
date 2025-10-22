from . import runner


def test_create_bbox_extracts():
    runner.run_snakemake_test(
        "create_bbox_extracts",
        (
            "results/json/djibouti-latest_extracts/slice-0.geojson",
            "results/json/djibouti-latest_extracts/slice-1.geojson",
            "results/json/djibouti-latest_extracts/slice-2.geojson",
            "results/json/djibouti-latest_extracts/slice-3.geojson",
        ),
    )
