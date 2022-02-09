# Take .geoparquet files and output a single, unified .geoparquet file
# TODO: recombining across hazards here means we haven't guaranteed that they
#   all use the same raster grid (because we checked consistently individually
#   for each hazard). We may need to tweak the flow to handle this.
rule join_data:
    input:
        expand(
            os.path.join(
                config['output_dir'],
                "splits",
                f"{dataset_slug}_filter-{filter_slug}_slice-{{i}}_hazard-{{hazard}}.geoparquet"
            ),
            i=range(config['slice_count']),
            hazard=config['hazard_datasets'].keys()
        ),
    output:
        # "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
        os.path.join(
            f"{config['output_dir']}",
            f"{dataset_slug}_filter-{filter_slug}.geoparquet"
        )
    script:
        "../scripts/join_data.py"

# This can be tested using
# snakemake --cores all join_data
