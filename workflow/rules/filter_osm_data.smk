# Take a .osm.pbf file and return a .osm.pbf file with a subset of the information
rule filter_osm_data:
    input:
        file="{OUTPUT_DIR}/input/{DATASET}.osm.pbf",
        filters=config["osmium_tags_filters_file"]
    output:
        "{OUTPUT_DIR}/input/{DATASET}_{FILTER_SLUG}.osm.pbf",
    shell:
        # filter for ways whose highway key has any of the values listed in filters
        "osmium tags-filter {input.file} w/highway=$(cat {input.filters}) -o {output}"

"""
Test with:
snakemake --cores all results/input/tanzania-mini_filter-highway-core.osm.pbf
"""
