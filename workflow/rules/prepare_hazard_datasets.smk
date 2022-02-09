# Create directories for each hazard in the config, and download its .txt file with resource locations to wget

rule prepare_hazard_datasets:
    output:
        expand(
            os.path.join(
                f"{config['output_dir']}",
                "input",
                f"hazard-{{hazard_name}}",
                f"file_list.txt"
            ),
            hazard_name=config['hazard_datasets'].keys()
        ),
    run:
        print(f"output={output}")
        print(f"wildcards={wildcards}")
        for hazard, txt_file in config['hazard_datasets'].items():
            if "_" in hazard or "/" in hazard:
                raise ValueError("Hazard names cannot contain _ or /")
            path = os.path.join(f"{config['output_dir']}", "input", f"hazard-{hazard}")
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            os.system(f"cd {path} && wget {txt_file} --output-document=file_list.txt")

rule test_prepare_hazard_datasets:
    input:
        expand(
            os.path.join(
                f"{config['output_dir']}",
                "input",
                f"hazard-{{hazard_name}}",
                f"{{hazard_name}}.txt"
            ),
            hazard_name=config['hazard_datasets'].keys()
        ),
