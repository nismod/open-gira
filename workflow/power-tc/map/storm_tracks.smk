rule plot_storm_tracks:
    """
    Draw map of storm tracks in set
    """
    input:
        tracks = "{OUTPUT_DIR}/storm_tracks/{STORM_SET}/tracks.geoparquet"
    params:
        storm_set_path = lambda wildcards: config["storm_sets"][wildcards.STORM_SET]
    output:
        png = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_tracks.png",
    run:
        import json
        import os

        import geopandas as gpd
        import matplotlib
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import numpy as np
        import pandas as pd
        import shapely

        from open_gira.wind import interpolate_track

        # MPL will segfault trying to plot as a snakemake subprocess otherwise
        matplotlib.use('Agg')

        tracks = gpd.read_parquet(input.tracks)

        with open(params.storm_set_path, "r") as fp:
            storm_set = json.load(fp)

        if len(storm_set) != 0:
            mask = tracks.track_id.isin(storm_set)
            data = tracks[mask]
        else:
            data = tracks

        # N.B. we use positive longitudes, so we can more easily center the map view
        # on the Pacific, with tracks spanning the antemeridian.

        data["longitude"] = data.geometry.x
        data["longitude"] = data["longitude"].apply(lambda l: l + 360 if l < 0 else l)
        data.geometry = gpd.points_from_xy(data.longitude, data.geometry.y)

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

        def rewrap_longitude(geometry: shapely.Geometry, meridian: float) -> shapely.Geometry:
            """
            Translate longitudes of geometries so they lie in the [meridian, meridian + 360] domain.

            Args:
                geometry: Geometric object to rewrap.
                meridian: Longitude to wrap at.

            Returns:
                Wrapped geometry.
            """
            return shapely.ops.transform(lambda x, y, z=None: (x + 360 if x < meridian else x, y), geometry)

        # Using longitude [0, 360], with wrapping on the Greenwich meridian
        # would break country boundary polygons. Instead, use longitude domain
        # [-20, 340], meaning that no land mass in the [-60, 60] latitude range
        # is split.
        world.geometry = world.geometry.apply(rewrap_longitude, args=(-20,))

#       interp = []
#       for track_id, groupby in data.groupby(data.track_id):
#           interp.append(interpolate_track(groupby, "1H"))
#       data = pd.concat(interp)

        f, ax = plt.subplots(figsize=(18,8))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1)

        ax = data.plot(
            column="max_wind_speed_ms",
            legend=True,
            legend_kwds={"label": "Wind speed [m s-1]"},
            vmin=20,
            vmax=80,
            cmap="plasma_r",
            alpha=0.2,
            s=1,
            ax=ax,
            cax=cax,
        )

        # draw borders
        world.boundary.plot(
            ax=ax,
            color="k",
            linewidth=1,
            alpha=0.2
        )

        ax.set_xlim(20, 340)
        ax.set_ylim(-60, 60)
        x_ticks = list(range(30, 360, 30))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{x:.0f}" if x <= 180 else f"{x-360:.0f}" for x in x_ticks])
        ax.grid(linestyle='--', alpha=0.3, linewidth=1)
        ax.set_xlabel("Longitude (degrees)")
        ax.set_ylabel("Latitude (degrees)")
        ax.set_title(f"{wildcards.STORM_SET} tropical cyclone tracks")

        f.savefig(output.png)

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/IBTrACS/storm_tracks.png
"""
