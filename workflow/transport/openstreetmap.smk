# Download the file specified in the config
import os
import re


rule download_osm:
    output:
        "{OUTPUT_DIR}/input/OSM/{DATASET}.osm.pbf",
    run:
        input_file = config["infrastructure_datasets"][wildcards.DATASET]
        if re.match("^https?://", input_file):
            os.system(f"wget {input_file} --output-document={output}")
        else:
            os.system("mkdir -p dirname {output}")
            os.system(f"cp {input_file} {output}")

        """
        Test with:
        snakemake --cores all results/input/OSM/tanzania-mini.osm.pbf
        """


# Take a .osm.pbf file and return a .osm.pbf file with a subset of the information
rule filter_osm_data:
    input:
        file="{OUTPUT_DIR}/input/OSM/{DATASET}.osm.pbf",
        # get the network type from the filter slug and use that to lookup the osmium filter expressions file
        filters=lambda wildcards: config["network_filters"][wildcards.FILTER_SLUG.replace("filter-", "")]
    output:
        "{OUTPUT_DIR}/input/OSM/{DATASET}_{FILTER_SLUG}.osm.pbf",
    shell:
        # filter for ways whose highway key has any of the values listed in filters
        "osmium tags-filter {input.file} --expressions {input.filters} -o {output}"


"""
Test with:
snakemake --cores all results/input/OSM/tanzania-mini_filter-road-primary.osm.pbf
"""
