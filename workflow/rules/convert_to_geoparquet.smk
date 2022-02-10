# Take .osm.pbf files and output .geoparquet files
rule convert_to_geoparquet:
    input:
        lambda wildcards: glob(
            f"{checkpoints.slice.get(**wildcards).output[0]}/{wildcards.SLICE_SLUG}.osm.pbf"
        ),
    output:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.geoparquet",
    script:
        "../scripts/osm_to_pq.py"

rule test_convert_to_geoparquet:
    input:
        expand(
            os.path.join(
                config['output_dir'],
                'geoparquet',
                f"{{dataset}}_slice-{{i}}_filter-{filter_slug}.geoparquet"
            ),
            dataset=config['infrastructure_datasets'].keys(),
            i=range(config['slice_count'])
        )