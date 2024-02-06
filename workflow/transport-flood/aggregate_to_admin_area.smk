rule aggregate_damages_to_admin_area:
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        slice_bounds = "{OUTPUT_DIR}/json/{DATASET}_extracts/{SLICE_SLUG}.geojson",
        data_to_aggregate = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{DIRECT_DAMAGE_TYPES}/{SLICE_SLUG}.geoparquet",
    params:
        columns_to_aggregate_regex = "^hazard-.+",
    output:
        aggregated_data = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{DIRECT_DAMAGE_TYPES}/{AGG_FUNC_SLUG}/{ADMIN_SLUG}/{SLICE_SLUG}.geoparquet",
    script:
        "./aggregate_to_admin_area.py"

"""
Test with:
snakemake -c1 results/direct_damages/europe-latest_filter-rail/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0/slice-0.geoparquet
"""


rule concat_and_sum_aggregated_damages:
    input:
        slices = lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "direct_damages",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "{DIRECT_DAMAGE_TYPES}",
                "{AGG_FUNC_SLUG}",
                "{ADMIN_SLUG}",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    params:
        columns_to_aggregate_regex = "^hazard-.+",
    output:
        joined = "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{DIRECT_DAMAGE_TYPES}/{AGG_FUNC_SLUG}/{ADMIN_SLUG}.geoparquet",
    script:
        "./concat_and_sum_slices.py"

"""
Test with:
snakemake --cores 1 results/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
"""
