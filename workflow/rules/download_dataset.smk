# Download the file specified in the config
rule download_dataset:
    output:
        os.path.join(f"{config['output_dir']}", "input", f"{dataset_slug}.osm.pbf"),
    shell:
        """
        wget {config[dataset]} --output-document={config[output_dir]}/input/{dataset_slug}.osm.pbf
        """

rule test_download_dataset:
    input:
        os.path.join(f"{config['output_dir']}", "input" ,f"{dataset_slug}.osm.pbf")