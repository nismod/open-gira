# Download the file specified in the config
rule download_dataset:
    output:
        "{OUTPUT_DIR}/input/{DATASET}.osm.pbf"
        # os.path.join(f"{config['output_dir']}", "input", f"{dataset_slug}.osm.pbf"),
    run:
        input_file = config['infrastructure_datasets'][wildcards.DATASET]
        os.system(f"wget {input_file} --output-document={output}")

"""
Test with:
snakemake --cores all results/input/tanzania-latest.osm.pbf
"""
