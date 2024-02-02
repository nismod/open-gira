"""Performs initial figure generation checks"""

import os

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
    models_all = snakemake.params["models_all"]  # type: ignore
except:
    raise RuntimeError("Please use snakemake to define inputs")


all_folders = [
    os.path.join(output_dir, f"power_output-{model}") for model in models_all
]

for folder in all_folders:  # check exist
    if not os.path.exists(folder):
        raise RuntimeWarning(
            f"\n----------------\nFolder {folder} does not exist. Can not proceed. Check that the workflow for {os.path.basename(folder)} has been completed and/or correctly named in {output_dir} directory.\n----------------\n"
        )

folder_figs = os.path.join(output_dir, "power_figures")
if not os.path.exists(folder_figs):
    os.makedirs(folder_figs)

txt_done = os.path.join(folder_figs, "initial_checks.txt")
with open(txt_done, "w") as f:
    f.writelines("Initial checks complete")
    print(f"written to {txt_done}")
