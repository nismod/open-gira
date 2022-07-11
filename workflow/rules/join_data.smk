# Take split .geoparquet files and output a single, unified .geoparquet file
# for each dataset-hazard combination


rule join_data:
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
        "../scripts/join_data.py"


"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet
"""
