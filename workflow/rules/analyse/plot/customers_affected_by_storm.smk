rule plot_top_storms_by_customers_affected:
    """
    Plot bar chart of customers affected by each storm, broken down by country
    """
    input:
        by_country = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption_by_country.nc",
        ibtracs = "{OUTPUT_DIR}/input/IBTrACS/processed/v4.geoparquet"
    output:
        pdf = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/customers_affected.pdf",
        pickle = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/customers_affected.pickle"
    run:
        import pickle

        import pandas as pd
        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import (FixedLocator, FixedFormatter)
        import xarray as xr

        MAX_STORMS = 20  # plot this many or fewer

        # MPL will segfault trying to plot as a snakemake subprocess otherwise
        matplotlib.use('Agg')

        # lookup from storm id to something human readable
        if "IBTrACS" in input.by_country:
            named_storms = True

            name_lookup = pd.read_parquet(input.ibtracs, columns=["year", "name", "track_id"]).set_index("track_id").drop_duplicates()
            name_lookup["name_year"] = name_lookup.name + ", " + name_lookup.year.apply(str)

        by_country: xr.Dataset = xr.open_dataset(input.by_country)

        data = {}
        min_threshold = min(by_country.threshold.values)
        for threshold in sorted(by_country.threshold.values):

            df = by_country.customers_affected.sel(dict(threshold=threshold)).to_pandas().T
            df = df.groupby(df.index).sum().sort_index()

            if threshold == min_threshold:
                # filter to most impactful storms, also drop all countries without impact
                # descending customers affected order
                storms_to_select = df.sum(axis=1).sort_values(ascending=False).index[:MAX_STORMS]
                df = df.loc[storms_to_select, :]
                countries_to_select = (df.sum(axis=0) != 0)
                df = df.loc[:, countries_to_select]

            # select all thresholds by same indicies as minimum threshold, order index by time
            data[threshold] = df.loc[storms_to_select, countries_to_select].sort_index()

        # make our own categorical colormap for the countries
        # allows us to persist country colours over several plots
        countries = list(data[min_threshold].columns)
        n_country = len(countries)
        colours = map(cm.cubehelix, np.linspace(0, 0.9, n_country))
        country_to_color = dict(zip(sorted(countries), colours))
        # table with RGBA rows and country columns
        colour_map = pd.DataFrame(country_to_color)

        # our x axes
        thresholds = list(data.keys())
        storm_ids = list(data[min_threshold].index.values)
        # the centre of each bar cluster on xaxis
        x = np.arange(len(storm_ids))

        width = 0.75
        l = range(len(thresholds))
        nudges = width / len(thresholds) * (l - np.median(l))

        bar_hatching = ["/" , r"\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ]

        f, ax = plt.subplots(figsize=(len(storm_ids) * 1 + 3, 10))

        # wind speed threshold loop
        minor_tick_locations = []
        for threshold, nudge in zip(thresholds, nudges):
            df: pd.DataFrame = data[threshold]
            top: pd.Series = df.cumsum(axis=1)
            bottom: pd.Series = top - df
            x_loc: np.ndarray = x + nudge

            # country loop
            for i, iso in enumerate(df.columns):
                colour = tuple(colour_map.loc[:, iso].values)
                bars = ax.bar(
                    x_loc,
                    df.loc[:, iso],
                    bottom=bottom.loc[:, iso],
                    log=False,
                    width=width / len(nudges),
                    label=iso,
                    color=colour,
                    alpha=0.5,
                    hatch=bar_hatching[i % len(bar_hatching)]
                )

            # store the minor tick locations
            minor_tick_locations.extend(x_loc)
#           ax.bar_label(
#               bars,
#               fmt='%.2G',
#               rotation=-90,
#               padding=3
#           )

        # axis labels and title
        subset_str = ""
        if len(by_country.event_id.values) > MAX_STORMS:
            subset_str = f" (top {MAX_STORMS} storms by customers affected)"
        ax.set_title(wildcards.STORM_SET + subset_str)
        ax.set_xlabel("Storm", labelpad=110)
        ax.set_ylabel("Customers affected", labelpad=40)
        ax.grid(which="both", alpha=0.2)

        # headroom for bar totals
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.1)

        # label the threshold speeds as minor ticks
        ax.set_xticks([], major=True)
        minor_tick_labels = list(thresholds) * len(storm_ids)
        minor_tick_locations = sorted(minor_tick_locations)
        assert len(minor_tick_labels) == len(minor_tick_locations)
        ax.xaxis.set_minor_formatter(FixedFormatter(minor_tick_labels))
        ax.xaxis.set_minor_locator(FixedLocator(minor_tick_locations))
        ax.xaxis.set_tick_params(which='minor', labelsize=8, rotation=-90)

        # draw the x-axis storm ids -- we're using the ticks for thresholds already
        for x_i, storm_id in zip(x, storm_ids):
            if named_storms:
                try:
                    name_str = name_lookup.loc[storm_id, 'name_year']
                except KeyError:
                    name_str = "Unknown"
                name_text = f"{name_str}\n{storm_id}"
            else:
                name_text = storm_id
            ax.text(
                x_i - width / len(thresholds),
                -0.17,  # below xaxis
                name_text,
                rotation=-90,
                size=8,
                horizontalalignment='left',
                verticalalignment='center',
                # transform: x axis in data units, y axis in fraction of plot
                transform=ax.get_xaxis_transform()
            )
        plt.subplots_adjust(bottom=0.22)
        ax.text(
            0.5,
            -0.085,
            "Wind speed threshold [m s-1]",
            horizontalalignment='center',
            # transform: x axis in data units, y axis in fraction of plot
            transform=ax.transAxes
        )

        # deduplicate and draw legend
        handles, labels = ax.get_legend_handles_labels()
        handle_list, label_list = [], []
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        ax.legend(
            handle_list,
            label_list,
            bbox_to_anchor=(1.02, 0.95),  # (x, y), anchor is top left of legend
            ncol=1,
            fancybox=True,
            shadow=True,
            handleheight=3,
            handlelength=4,
        )

        with open(output.pickle, "wb") as fp:
            pickle.dump(f, fp)

        f.savefig(output.pdf)

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS-maria_2017/customers_affected.pdf
"""
