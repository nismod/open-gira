"""
Map wind fields
"""


rule map_wind_field:
    """
    Draw a map of the merged wind field for a given storm.
    """
    input:
        wind_field = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/wind_field.nc",
    output:
        plot = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/wind_field.png",
    run:
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        import geopandas as gpd
        import xarray as xr

        from open_gira.plot.utils import figure_size

        # MPL will segfault trying to plot as a snakemake subprocess otherwise
        matplotlib.use('Agg')

        ds = xr.open_dataset(input.wind_field)

        if not ds.max_wind_speed.any():
            # no data in netCDF
            f, ax = plt.subplots()
            ax.set_title(wildcards.STORM_ID)
            f.savefig(output.plot)

        else:
            f, ax = plt.subplots(
                figsize=figure_size(
                    ds.longitude.min(),
                    ds.latitude.min(),
                    ds.longitude.max(),
                    ds.latitude.max()
                )
            )
            plt.subplots_adjust(left=0.1, right=0.8)
            cax = f.add_axes(
                # x_min, y_min, x_delta, y_delta
                [0.84, 0.11, 0.04, 0.77]
            )
            vmin = 21
            vmax = 45
            levels = int((vmax - vmin) / 3 + 1)
            xr.plot.pcolormesh(
                ds.max_wind_speed,
                x="longitude",
                y="latitude",
                extend="max",
                cmap="turbo",
                vmin=vmin,
                vmax=vmax,
                levels=levels,
                ax=ax,
                cbar_ax=cax,
                alpha=0.7
            )
            cax.set_ylabel("Maximum wind speed [ms-1]", labelpad=10)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_title(wildcards.STORM_ID)

            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))

            ax.grid()

            borders = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            borders.plot(ax=ax, facecolor="none", edgecolor="grey", alpha=0.8)

            f.savefig(output.plot)

"""
To test:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017260N12310/wind_field.png
"""


def wind_field_map_paths_for_storm_set(wildcards):
    """
    Generate list of paths for wind fields in storm set
    """
    import json
    import os

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    with open(json_file, "r") as fp:
        country_set_by_storm = json.load(fp)

    paths = []
    for storm_id in country_set_by_storm.keys():
        paths.extend(
            expand(
                "results/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/wind_field.png",
                STORM_SET=wildcards.STORM_SET,  # str
                STORM_ID=storm_id  # str
            )
        )

    return paths


rule map_storm_set:
    """
    Draw maps of all wind fields for storm set and symlink these files into one directory.
    """
    input:
        wind_field_maps = wind_field_map_paths_for_storm_set
    output:
        map_dir = directory("{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/wind_field_maps")
    shell:
        """
        mkdir -p {output.map_dir}
        for FILE in {input.wind_field_maps}; do
            STORM_ID=$(echo $FILE | cut -d/ -f6)
            ln -s $(pwd)/$FILE $(pwd)/{output.map_dir}/${{STORM_ID}}.png
        done
        """

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/wind_field_maps
"""
