"""Plots the cdfs for EAD and EACA

"""
import os

cdf_EACA = "empirical_effective population affected_plotting_data"
cdf_EAD = "empirical_reconstruction cost_plotting_data"
in_cdf_EAD = [
    os.path.join(
        config["output_dir"],
        f"power_output-{model}",
        "statistics",
        "empirical",
        "empirical_plotting_data",
        f"{cdf_EAD}.csv",
    )
    for model in models_all
]
in_cdf_EACA = [
    os.path.join(
        config["output_dir"],
        f"power_output-{model}",
        "statistics",
        "empirical",
        "empirical_plotting_data",
        f"{cdf_EACA}.csv",
    )
    for model in models_all
]
out_cdf_EAD = os.path.join(
    config["output_dir"], "power_figures", f"empirical_overlay_{cdf_EAD}.png"
)
out_cdf_EACA = os.path.join(
    config["output_dir"], "power_figures", f"empirical_overlay_{cdf_EACA}.png"
)


rule fig_cdfs_EAD:
    conda: "../../../../environment.yml"
    input:
        in_cdf_EAD,
    params:
        output_dir=config["output_dir"],
        EACA=False,
        central_threshold=config["central_threshold"],
        minimum_threshold=config["minimum_threshold"],
        maximum_threshold=config["maximum_threshold"],
    output:
        out_cdf_EAD,
    script:
        os.path.join(
            "..", "..", "..", "scripts", "analyse", "figures", "plot_together.py"
        )


rule fig_cdfs_EACA:
    conda: "../../../../environment.yml"
    input:
        in_cdf_EACA,
    params:
        output_dir=config["output_dir"],
        EACA=True,  # EAD
        central_threshold=config["central_threshold"],
        minimum_threshold=config["minimum_threshold"],
        maximum_threshold=config["maximum_threshold"],
    output:
        out_cdf_EACA,
    script:
        os.path.join(
            "..", "..", "..", "scripts", "analyse", "figures", "plot_together.py"
        )
