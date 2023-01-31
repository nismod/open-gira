rule plot_storm_tracks:
    """
    Draw map of storm tracks in set
    """
    input:
        storm_set = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_set.json",
        ibtracs = "{OUTPUT_DIR}/input/IBTrACS/processed/v4.geoparquet"
    output:
        pdf = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_tracks.pdf",
        pickle = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_tracks.pickle"
    run:
        import json
        import os
        import pickle

        import numpy as np
        import pandas as pd
        import geopandas as gpd
        import matplotlib
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        from open_gira.wind import interpolate_track

        # MPL will segfault trying to plot as a snakemake subprocess otherwise
        matplotlib.use('Agg')

        ibtracs = gpd.read_parquet(input.ibtracs)
        ibtracs["name_year"] = ibtracs["name"] + " " + ibtracs["year"].apply(str)

        with open(input.storm_set, "r") as fp:
            storm_set = json.load(fp)

        mask = ibtracs.track_id.isin(storm_set)
        data = ibtracs[mask]

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

        interp = []
        for track_id, groupby in data.groupby(data.track_id):
            interp.append(interpolate_track(groupby, "1H"))
        interp = pd.concat(interp)

        f, ax = plt.subplots(figsize=(20,10))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cax.set_title("Wind speed [m s-1]")

        ax = interp.plot(
            column="max_wind_speed_ms",
            figsize=(18, 12),
            legend=True,
            s=4,
            ax=ax,
            cax=cax,
        )

        # draw borders
        world.boundary.plot(
            ax=ax,
            color="k",
            alpha=0.2
        )

        ax.set_xlim(-180, 180)
        ax.set_ylim(-80, 80)

        ax.set_title(f"{wildcards.STORM_SET} tracks")

        # label storms
        for track_id, track in interp.groupby(interp.track_id):
            ax.annotate(
                text=track.name_year.iloc[0],
                fontsize=10,
                xy=track.geometry.iloc[0].coords[0],
                ha='left',
                rotation=0
            )

        with open(output.pickle, "wb") as fp:
            pickle.dump(f, fp)

        f.savefig(output.pdf)
