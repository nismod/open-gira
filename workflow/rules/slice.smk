# Use osmium to cut a .osm.pbf file into several .osm.pbf files as defined by bounding boxes
rule slice:
    input:
        # data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        data=os.path.join(f"{config['output_dir']}",f"{config['dataset']}_filter-{filter_slug}.osm.pbf"),
        extracts_config=os.path.join(
            f"{config['output_dir']}",
            "json",
            f"{config['dataset']}_filter-{filter_slug}-extracts.geojson"
        ),
    output:
        "{OUTPUT_DIR}/slices/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.osm.pbf",
    shell:
        """
        osmium extract --overwrite --no-progress --config {input.extracts_config} {input.data}
        """

rule test_slice:
    input:
        expand(
            os.path.join(
                f"{config['output_dir']}",
                "slices",
                f"{config['dataset']}_slice-{{i}}_filter-{filter_slug}.osm.pbf"
            ),
            i=range(config['slice_count'])
        ),