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
