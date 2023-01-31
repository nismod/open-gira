rule plot_target_supply_factor_distributions:
    """
    Plot distributions of supply_factor, from the pool of targets processed for storm_set
    """
    input:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure_by_target.nc",
    output:
        pdf = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/supply_factor_distribution.pdf",
    run:
        import matplotlib
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        import numpy as np
        import xarray as xr

        # MPL will segfault trying to plot as a snakemake subprocess otherwise
        matplotlib.use('Agg')

        ds = xr.open_dataset(input.by_target)

        # preprocess the data for plotting
        data = {}
        for threshold in ds.threshold.values:
            # flatten
            supply_factor = ds.supply_factor.sel(dict(threshold=threshold)).values.ravel()
            # drop >= 1, and any NaN values
            supply_factor = supply_factor[~np.isnan(supply_factor)]
            supply_factor = supply_factor[supply_factor < 1]

            data[threshold] = supply_factor

        n_threshold = len(ds.threshold.values)

        f, axes = plt.subplots(nrows=1, ncols=n_threshold, sharex=True, figsize=(2 + 2 * n_threshold, 4))

        colours = list(mcolors.BASE_COLORS.values())
        min_y_all = 1
        max_y_all = 0
        for i, (threshold, supply_factor) in enumerate(data.items()):
            axes[i].hist(
                supply_factor,
                bins=25,
                label=threshold,
                alpha=0.5,
                density=True,
                color=colours[i % len(colours)]
            )

            min_y, max_y = axes[i].get_ylim()
            min_y_all = min(min_y, min_y_all)
            max_y_all = max(max_y, max_y_all)

            axes[i].set_xlim(0, 1)
            axes[i].grid(which="both", alpha=0.2)
            if i != 0:
                axes[i].set_yticklabels([])

        for ax in axes:
            if min_y_all != 0:
                axes[i].set_yscale("log")
            ax.set_ylim(0.5 * min_y_all, 1.1 * max_y_all)

        f.suptitle("Distribution of supply factor < 1")
        f.supxlabel("Supply factor")
        f.supylabel("Density")

        lines_labels = [ax.get_legend_handles_labels() for ax in axes]
        lines, labels = [sum(ll, []) for ll in zip(*lines_labels)]
        f.legend(
            lines,
            labels,
            title="Damage threshold [m s-1]"
        )

        plt.subplots_adjust(bottom=0.15)

        f.savefig(output.pdf)
