# Take a .osm.pbf file and return a .osm.pbf file with a subset of the information
rule filter_osm_data:
    input:
        config["osmium_tags_filters_file"],
        "{OUTPUT_DIR}/input/{DATASET}.osm.pbf"
    output:
        "{OUTPUT_DIR}/input/{DATASET}_{FILTER_SLUG}.osm.pbf",
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"

"""
Test with:
snakemake --cores all results/input/tanzania-latest_filter-highway-core.osm.pbf
"""
