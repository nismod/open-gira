# Take .geoparquet files and output a single, unified .geoparquet file
rule join_data:
    input:
        expand(
            os.path.join(
                f"{config['output_dir']}",
                "splits",
                f"{config['dataset']}_filter-{filter_slug}_slice-{{i}}_hazard-{hazard_slug}.geoparquet"
            ),
            i=range(config['slice_count'])
        ),
    output:
        f"{config['output_dir']}/{config['dataset']}_filter-{filter_slug}_hazard-{hazard_slug}.geoparquet",
    script:
        "../scripts/join_data.py"

rule test_join_data:
    input:
        os.path.join(
            config['output_dir'],
            f"{config['dataset']}_filter-{filter_slug}_hazard-{hazard_slug}.geoparquet"
        ),