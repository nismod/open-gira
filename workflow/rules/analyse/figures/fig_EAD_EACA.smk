"""Writes EAD and EACA metrics to file in power_figures

"""
import os

out_EAD_EACA = [
    os.path.join(config["output_dir"], "power_figures", f"{metric}_total.txt")
    for metric in ["EAD", "EACA"]
]


rule fig_EAD_EACA:
    input:
        [
            os.path.join(
                config["output_dir"],
                f"power_output-{model}",
                "statistics",
                "empirical",
                "empirical_plotting_data",
                "empirical_effective population affected_plotting_data.csv",
            )
            for model in models_all
        ],
        [
            os.path.join(
                config["output_dir"],
                f"power_output-{model}",
                "statistics",
                "empirical",
                "empirical_plotting_data",
                "empirical_reconstruction cost_plotting_data.csv",
            )
            for model in models_all
        ],
    params:
        output_dir=config["output_dir"],
        models_future=models_future,
    output:
        out_EAD_EACA,
    script:
        os.path.join(
            "..", "..", "..", "scripts", "analyse", "figures", "EAD_EACA-finder.py"
        )
