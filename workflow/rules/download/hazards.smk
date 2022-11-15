# Download hazard data from the file_list.txt files in each hazard directory
import os
import re

# This is a checkpoint because it outputs to a directory we want to ensure is up to date in later rules.
# Specifically, trim_hazard_data.smk looks for all *.tif files in the output directory.
from tempfile import TemporaryDirectory


checkpoint download_hazard_datasets:
    output:
        raw_folder=directory("{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw"),
    run:
        with TemporaryDirectory() as tmpdir:
            input_file_key = re.sub("^hazard-", "", wildcards.HAZARD_SLUG)
            input_file = config["hazard_datasets"][input_file_key]
            # Download remote input file if necessary
            if re.match("https?://", input_file):
                os.system(
                    f"wget {input_file} --output-document={tmpdir}/hazard_sources.txt"
                )
                input_file = f"{tmpdir}/hazard_sources.txt"

            # Split input_file into local and remote resources
            local_files = []
            remote_files = []
            with open(input_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if re.match("^https?://", line):
                        remote_files.append(f"{line}\n")
                    else:
                        local_files.append(line)

            print(f"local_files={local_files}")
            print(f"remote_files={remote_files}")
            target_dir = f"{wildcards.OUTPUT_DIR}/input/{wildcards.HAZARD_SLUG}/raw"
            os.system(f"mkdir --parents {target_dir}")

            for f in local_files:
                os.system(f"cp {f} {target_dir}/{os.path.basename(f)}")

            if len(remote_files):
                with open(f"{tmpdir}/input.txt", "w") as sources:
                    sources.writelines("\n".join(remote_files))
                # use xargs to run multiple wget processes, one per hazard URL
                # -n 1 to expect only one argument per process
                # -P {workflow.cores} to limit the total concurrent connections to --cores value
                os.system(f"cd {target_dir} && xargs -n 1 -P {workflow.cores} wget --no-clobber < {tmpdir}/input.txt")

        """
        Test with:
        snakemake --cores all results/input/hazard-aqueduct-river/raw
        """
