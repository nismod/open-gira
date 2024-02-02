rule join_data:
    """
    Take split .geoparquet files and output a single, unified .geoparquet file
    for each dataset-hazard combination
    """
    input:
        lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "splits",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
    script:
        "./join_data.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-road_hazard-aqueduct-river.geoparquet
"""

# Take split .geoparquet files and output a single, unified .geoparquet file
# for each dataset-hazard combination


rule join_exposure:
    input:
        lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "splits",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
    script:
        "./join_data.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet
"""


rule join_direct_damages:
    input:
        slices = lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "direct_damages",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "{DIRECT_DAMAGE_TYPES}",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        joined = "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_{DIRECT_DAMAGE_TYPES}.geoparquet",
    script:
        "./join_data.py"

"""
Test with:
snakemake --cores all results/egypt-latest_filter-road/hazard-aqueduct-river/damage_fraction_per_RP.geoparquet
"""
