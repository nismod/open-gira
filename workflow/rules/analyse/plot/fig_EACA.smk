"""
Aggregates EACA and finds difference to future. Plots current and difference
"""
import os

increased_severity_sort_bool = str(config["increased_severity_sort"])[0]
assert increased_severity_sort_bool in ["T", "F"]

metric_EACA = "effective_population_anually-expected"
metric_EACA_perc = "perc_effective_population_anually-expected"
merge_key_EACA = "code"

out_agg_EACA_file = os.path.join(
    config["output_dir"],
    "power_figures",
    "intermediate_files",
    f"mean_{merge_key_EACA}_{metric_EACA}.gpkg",
)
out_agg_EACA_plot_perc = os.path.join(
    config["output_dir"], "power_figures", "EACA_perc_diff_aggregate.png"
)
out_agg_EACA_plot = os.path.join(
    config["output_dir"], "power_figures", "EACA_current_aggregate.png"
)

out_diff_EACA_file = os.path.join(
    config["output_dir"],
    "power_figures",
    "intermediate_files",
    f"difference_{merge_key_EACA}_{metric_EACA}.gpkg",
)


rule fig_aggregate_EACA:
    input:
        in_agg_EACA=[
            os.path.join(
                config["output_dir"],
                f"power_output-{model}",
                "statistics",
                "aggregate",
                f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent_aggregated_region.gpkg",
            )
            for model in models_future
        ],
    params:
        output_dir=config["output_dir"],
        metric=metric_EACA,
        merge_key=merge_key_EACA,
    output:
        out_agg_EACA=out_agg_EACA_file,  #
    script:
        os.path.join("..", "..", "..", "scripts", "analyse", "figures", "mean_agg.py")


rule fig_diff_EACA:
    input:
        [
            rules.fig_aggregate_EACA.output.out_agg_EACA,
            os.path.join(
                config["output_dir"],
                f"power_output-constant",
                "statistics",
                "aggregate",
                f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent_aggregated_region.gpkg",
            ),
        ],  # first file is future
    params:
        output_dir=config["output_dir"],
        metric=metric_EACA,
        merge_key=merge_key_EACA,
    output:
        out_diff_EACA=out_diff_EACA_file,
    script:
        os.path.join("..", "..", "..", "scripts", "analyse", "figures", "diff_agg.py")


rule fig_plot_current_EACA:
    """Plots current"""
    input:
        rules.fig_aggregate_EACA.input.in_agg_EACA[0],  # constant
    params:
        output_dir=config["output_dir"],
        metric=metric_EACA,
        vmax=50000,
        vmin=0,
        cmap="Reds",
        linewidth=None,
        legend_name="EACA (constant climate) [people]",
    output:
        out_agg_EACA_plot,
    script:
        os.path.join("..", "..", "..", "scripts", "analyse", "figures", "plotter.py")


rule fig_plot_diff_EACA:
    """Plots difference"""
    input:
        rules.fig_diff_EACA.output.out_diff_EACA,
    params:
        output_dir=config["output_dir"],
        metric=metric_EACA_perc,
        vmax=30,
        vmin=-30,
        cmap="RdBu_r",
        linewidth=None,
        legend_name="EACA change [% of constant climate]",
    output:
        out_agg_EACA_plot_perc,
    script:
        os.path.join("..", "..", "..", "scripts", "analyse", "figures", "plotter.py")
