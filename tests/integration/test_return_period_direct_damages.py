from . import runner


def test_return_period_direct_damages():
    runner.run_snakemake_test(
        "return_period_direct_damages",
        (
            "results/direct_damages/djibouti-latest_filter-road/hazard-aqueduct-river/fraction_per_RP/slice-0.geoparquet",
            "results/direct_damages/djibouti-latest_filter-road/hazard-aqueduct-river/cost_per_RP/slice-0.geoparquet",
            "results/direct_damages/djibouti-latest_filter-road/hazard-aqueduct-river/EAD/slice-0.geoparquet",
            "results/direct_damages/djibouti-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/slice-0.geoparquet",
        )
    )
