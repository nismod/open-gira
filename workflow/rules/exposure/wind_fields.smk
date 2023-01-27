"""
Estimate max wind speed at infrastructure asset locations per event
"""


rule create_wind_grid:
    """
    Create an empty TIFF file for a given box specifying the spatial grid to
    evaluate wind speed on
    """
    input:
        network_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json",
    output:
        wind_grid="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/wind_grid.tiff",
    run:
        import os
        import json
        import sys

        import shapely
        import numpy as np

        def harmonise_grid(minimum: float, maximum: float, cell_length: float) -> tuple[int, float, float]:
            """
            Grow grid dimensions to encompass whole number of `cell_length`

            Args:
                minimum: Minimum dimension value
                maximum: Maximum dimension value
                cell_length: Length of cell side

            Returns:
                Number of cells
                Adjusted minimum
                Adjusted maximum
            """

            assert maximum > minimum

            span: float = maximum - minimum
            n_cells: int = int(np.ceil(span / cell_length))
            delta: float = n_cells * cell_length - span
            buffer: float = delta / 2

            return n_cells, minimum - buffer, maximum + buffer

        # read hull shape from disk
        with open(input.network_hull, "r") as fp:
            data = json.load(fp)

        try:
            shape_dict, = data["features"]
        except ValueError:
            # no grid, touch empty file and exit
            os.system(f"touch {output.wind_grid}")
            return

        hull = shapely.geometry.shape(shape_dict["geometry"])

        minx, miny, maxx, maxy = hull.bounds
        cell_length = config["wind_deg"]  # cell side length in decimal degrees

        # determine grid bounding box to fit an integer number of grid cells in each dimension
        i, minx, maxx = harmonise_grid(minx, maxx, cell_length)
        j, miny, maxy = harmonise_grid(miny, maxy, cell_length)

        # create grid as TIFF and save to disk
        os.makedirs(os.path.dirname(output.wind_grid), exist_ok=True)
        command = f"gdal_create -outsize {i} {j} -a_srs EPSG:4326 -a_ullr {minx} {miny} {maxx} {maxy} {output.wind_grid}"
        os.system(command)

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/storms/wind_grid.tiff
"""


def storm_dataset(wildcards) -> str:
    """Return path to storm tracks"""
    # parent dataset of storm set e.g. IBTrACS for IBTrACS_maria-2017
    storm_dataset = wildcards.STORM_SET.split("_")[0]
    return f"{wildcards.OUTPUT_DIR}/power/by_country/{wildcards.COUNTRY_ISO_A3}/storms/{storm_dataset}/tracks.geoparquet"


def read_storm_set(wildcards) -> list[str]:
    """
    Read file containing a list of storm IDs that constitute the storm set.

    N.B. An empty file (and set) will be interpreted as process all storms in
    the dataset.
    """
    storm_set_path = config["storm_sets"][wildcards.STORM_SET]
    with open(storm_set_path, "r") as fp:
        storm_set = json.load(fp)
    return storm_set


rule estimate_wind_fields:
    """
    Find maximum windspeeds for each storm for each grid cell.

    Optionally plot wind fields and save to disk
    """
    conda: "../../../environment.yml"
    input:
        storm_file=storm_dataset,
        wind_grid="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/wind_grid.tiff",
    params:
        storm_set=read_storm_set
    threads:
        config["processes_per_parallel_job"]
    output:
        # enable or disable plotting in the config file
        plot_dir=directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/plots/"),
        wind_speeds="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/max_wind_field.nc",
    script:
        "../../scripts/intersect/estimate_wind_fields.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IBTrACS/max_wind_field.nc
"""
