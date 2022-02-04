# Take .geoparquet files and output a single, unified .geoparquet file
rule join_data:
    input:
        expand(
            os.path.join(
                f"{config['output_dir']}",
                "splits",
                f"{config['dataset']}-slice{{i}}.highway-core_{hazard_slug}_splits.geoparquet"
            ),
            i=range(config['slice_count'])
        ),
    output:
        f"{config['output_dir']}/{config['dataset']}.highway-core_{hazard_slug}_splits.geoparquet",
    script:
        "../scripts/join_data.py"
