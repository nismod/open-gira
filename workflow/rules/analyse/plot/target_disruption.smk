rule plot_target_supply_factor_distributions:
    """
    Plot distributions of supply_factor, from the pool of targets processed for storm_set
    """
    input:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption_by_target.nc",
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

        n_threshold = len(ds.threshold.values)

        f, axes = plt.subplots(nrows=1, ncols=n_threshold, sharex=True, figsize=(2 + 2 * n_threshold, 4))

        colours: list[tuple[float]] = [
            matplotlib.colormaps["cubehelix"](i) for i in np.linspace(0, 0.9, n_threshold)
        ]

        max_y_all = 0
        y_cutoff = 1
        for i, threshold in enumerate(ds.threshold.values):

            # pool values across storms and targets
            supply = ds.supply_factor.sel(dict(threshold=threshold)).values.ravel()
            # drop NaN values
            supply = supply[~np.isnan(supply)]
            # bin the data
            freq, edges = np.histogram(supply, bins=np.linspace(0, 1, 11))
            # set 0 values to 1 (permits log)
            freq[freq <= 0] = y_cutoff
            width = edges[1] - edges[0]
            axes[i].bar(
                edges[:-1] + width / 2,
                freq,
                width=width,
                label=threshold,
                alpha=0.8,
                color=colours[i]
            )

            _, max_y = axes[i].get_ylim()
            max_y_all = max(max_y, max_y_all)

            axes[i].set_yscale("log")
            axes[i].set_xlim(0, 1)
            axes[i].grid(which="both", alpha=0.2)
            axes[i].set_title(threshold)
            if i != 0:
                axes[i].set_yticklabels([])

        for ax in axes:
            ax.set_ylim(y_cutoff, 2 * max_y_all)

        f.suptitle("Supply factor by damage threshold [ms-1]")
        f.supxlabel("Supply factor")
        f.supylabel("Frequency")

        plt.subplots_adjust(bottom=0.15, top=0.85, left=0.075, right=0.95)

        f.savefig(output.pdf)
