# Take split .geoparquet files and output a single, unified .geoparquet file
# for each dataset-hazard combination

rule join_data:
    input:
        lambda wildcards: expand(
            os.path.join("{OUTPUT_DIR}", "splits", "{DATASET}_{FILTER_SLUG}_slice-{i}_{HAZARD_SLUG}.geoparquet"),
            **wildcards,
            i=range(config['slice_count'])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
        # os.path.join(
        #     f"{config['output_dir']}",
        #     f"{dataset_slug}_filter-{filter_slug}.geoparquet"
        # )
    script:
        "../scripts/join_data.py"

"""
Test with:
snakemake --cores all results/tanzania-latest_filter-highway-core_hazard-aqueduct-coast.geoparquet
"""
