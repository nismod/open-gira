"""
Estimate max wind speed at infrastructure asset locations per event
"""

from open_gira.io import cached_json_file_read


rule create_wind_grid:
    """
    Create an empty TIFF file for a given box specifying the spatial grid to
    evaluate wind speed on
    """
    input:
        network_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json",
    params:
        # include wind_grid_resolution_deg as a param to trigger re-runs on change
        grid_resolution=config["wind_grid_resolution_deg"]
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

        # expand grid by buffer_deg, gives us a view of storm intensity over coastal zone
        # model grid slightly smaller than track search radius, should mean storm
        # eyes cross boundary rather than appearing unannounced
        buffer_deg = config["max_track_search_radius_deg"] - 0.5
        minx, miny, maxx, maxy = hull.bounds
        minx -= buffer_deg
        miny -= buffer_deg
        maxx += buffer_deg
        maxy += buffer_deg

        # cell side length in decimal degrees
        cell_length = config["wind_grid_resolution_deg"]

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


rule cut_land_cover_raster:
    """
    Crop global land cover raster to country.
    """
    input:
        global_raster = "{OUTPUT_DIR}/input/land_cover/glob_cover_2009/GLOBCOVER_L4_200901_200912_V2.3.tif",
        wind_grid = rules.create_wind_grid.output.wind_grid
    output:
        local_raster = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/land_cover.tiff",
    shell:
        """
        # the wind grid file may be empty (we couldn't make a network, grid hull file and therefore wind grid)
        # in that case, gdalinfo will fail, so check the exit status and branch if necessary
        # must disable snakemake's bash strict mode to permit a non-zero exit code
        set +e
        DATA=$(gdalinfo -json {input.wind_grid})
        if [ "$?" -ne 0 ]
        then
            echo "Found empty wind grid file, writing empty land cover raster"
            touch {output.local_raster}
            exit 0
        fi
        # re-enable strict mode
        set -e

        # N.B. upper and lower are the reverse of what you'd expect
        # i.e. the minimum latitude is at the upper corners
        # and the maximum latitude is at the lower corners
        # longitude as expected
        MINX=$(echo $DATA | jq .cornerCoordinates.lowerLeft[0])
        MINY=$(echo $DATA | jq .cornerCoordinates.upperRight[1])
        MAXY=$(echo $DATA | jq .cornerCoordinates.lowerLeft[1])
        MAXX=$(echo $DATA | jq .cornerCoordinates.upperRight[0])

        echo "Cropping to following bounding box..."
        echo $MINX $MINY $MAXX $MAXY
        gdalwarp -te $MINX $MINY $MAXX $MAXY {input.global_raster} {output.local_raster}
        """

"""
Test with:
snakemake -c1 results/power/by_country/CUB/storms/land_cover.tiff
"""


rule create_surface_roughness_raster:
    """
    1) Load the high resolution (300m) land cover classification
    2) Lookup surface roughness values
    3) Downsample to native wind grid
    """
    input:
        land_cover = rules.cut_land_cover_raster.output.local_raster,
        wind_grid = rules.create_wind_grid.output.wind_grid,
        land_cover_roughness_mapping = "config/land_use_to_surface_roughness.csv"
    output:
        surface_roughness = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/surface_roughness.tiff"
    script:
        "./surface_roughness.py"

"""
Test with:
snakemake -c1 results/power/by_country/CUB/surface_roughness.tiff
"""


rule create_downscaling_factors:
    """
    Use surface roughness values calculate factors for scaling from gradient-level to surface-level winds.
    """
    input:
        surface_roughness=rules.create_surface_roughness_raster.output.surface_roughness,
    output:
        downscale_factors="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/downscale_factors.npy",
        downscale_factors_plot="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/downscale_factors.png",
    script:
        "./wind_downscaling_factors.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/downscale_factors.npy
"""


def storm_tracks_by_country(wildcards) -> str:
    """Return path to storm tracks"""
    # parent dataset of storm set e.g. IBTrACS for IBTrACS_maria-2017
    storm_dataset = wildcards.STORM_SET.split("_")[0]
    return f"{wildcards.OUTPUT_DIR}/power/by_country/{wildcards.COUNTRY_ISO_A3}/storms/{storm_dataset}/{wildcards.SAMPLE}/tracks.geoparquet"


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
    input:
        # `threads_for_country` will fail unless this CSV is present when resources are set
        country_target_count=country_target_count_path,
        storm_file=storm_tracks_by_country,
        wind_grid="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/wind_grid.tiff",
        downscaling_factors=rules.create_downscaling_factors.output.downscale_factors,
    params:
        storm_set=read_storm_set,
        # include failure_thresholds as a param (despite not using elsewhere in the
        # rule) to trigger re-runs on change to this configuration option
        failure_thresholds=config["transmission_windspeed_failure"]
    threads: threads_for_country
    resources:
        # 2GB RAM per CPU
        mem_mb = lambda wildcards: threads_for_country(wildcards) * 2_048
    output:
        # enable or disable plotting in the config file
        plot_dir=directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/{SAMPLE}/plots/"),
        wind_speeds="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/{SAMPLE}/max_wind_field.nc",
    script:
        "./estimate_wind_fields.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IBTrACS/0/max_wind_field.nc
"""


def wind_field_paths_all_samples(wildcards) -> list[str]:
    """
    Return a list of paths for every sample
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/{SAMPLE}/max_wind_field.nc",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name])
    )


rule concat_wind_fields_over_sample:
    """
    Take wind fields expanded over sample and stack them into a single netCDF.
    """
    input:
        sample_paths=wind_field_paths_all_samples
    output:
        concat="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/max_wind_field.nc",
    script:
        "./concat_wind_over_sample.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/STORM-constant/max_wind_field.nc
"""


def wind_fields_by_country_for_storm(wildcards):
    """
    Given a STORM_ID and STORM_SET as a wildcard, lookup the countries that storm
    affects and return paths to their wind field files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    country_set_by_storm = cached_json_file_read(json_file)

    return expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/{SAMPLE}/max_wind_field.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        SAMPLE=wildcards.SAMPLE  # str
    )


rule merge_wind_fields_of_storm:
    """
    Merge wind fields generated for each country for a given storm.
    """
    input:
        wind_fields = wind_fields_by_country_for_storm
    output:
        merged = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/{SAMPLE}/by_storm/{STORM_ID}/wind_field.nc",
    run:
        import logging
        import os

        import geopandas as gpd
        import pandas as pd
        import numpy as np
        import xarray as xr

        from open_gira.wind import empty_wind_da
        from open_gira.binning import grid_point_data

        # filter out really low wind speeds to constrain the map extent
        MIN_WIND_SPEED = 10

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Reading wind fields for each country")
        rasters = []
        for path in input.wind_fields:
            try:
                ds = xr.open_dataset(path).sel(event_id=wildcards.STORM_ID)
                rasters.append(ds.to_dataframe().reset_index()[["latitude", "longitude", "max_wind_speed"]])
            except KeyError:
                pass

        try:
            data = pd.concat(rasters)
        except ValueError:
            # no data in rasters
            data = pd.DataFrame({"max_wind_speed": []})

        logging.info(f"Filtering out areas with wind speed < {MIN_WIND_SPEED} ms-1")
        data = data[data.max_wind_speed > MIN_WIND_SPEED]

        if len(data) != 0:
            # transform to point data
            delta = float(config["wind_grid_resolution_deg"])
            logging.info(f"Regridding on harmonised {delta} degree grid.")
            point_speeds = gpd.GeoDataFrame(
                {
                    "max_wind_speed": data.max_wind_speed,
                    "geometry": gpd.points_from_xy(data.longitude, data.latitude)
                }
            )
            # rebin on new grid
            rebinned: gpd.GeoDataFrame = grid_point_data(point_speeds, "max_wind_speed", "max", delta)

            # transform to xarray
            raster = rebinned.copy()

            # round coordinates to 6 decimal places to remove coordinate gaps arising from floating point (im)precision
            # this is necessary to remove spurious extra coordinates from appearing in output netCDF
            raster["longitude"] = raster.geometry.centroid.x.round(6)
            raster["latitude"] = raster.geometry.centroid.y.round(6)

            raster = raster.drop(columns=["geometry"])
            raster: gpd.GeoDataFrame = raster.set_index(["longitude", "latitude"])
            rebinned_raster: xr.Dataset = raster.to_xarray()

            rebinned_raster.to_netcdf(output.merged)

        else:
            logging.info(f"No data in rasters, writing empty file.")
            empty_wind_da().to_netcdf(output.merged)

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS/0/by_storm/2017260N12310/wind_field.nc
"""


def merged_wind_fields_for_all_storms_in_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the storms in the set.

    Return a list of the merged wind_field.nc file paths for every storm in the set.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    country_set_by_storm = cached_json_file_read(json_file)

    storms = list(country_set_by_storm.keys())

    return expand(
        "results/power/by_storm_set/{STORM_SET}/{SAMPLE}/by_storm/{STORM_ID}/wind_field.nc",
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=storms,  # list of str
        SAMPLE=wildcards.SAMPLE  # str
    )


rule merged_wind_fields_for_storm_set:
    """
    A target rule to generate the wind field netCDFs (across multiple
    countries) for each storm.
    """
    input:
        wind_fields = merged_wind_fields_for_all_storms_in_storm_set
    output:
        completion_flag = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/{SAMPLE}/wind_fields.txt"
    shell:
        """
        # one output file per line
        echo {input.wind_fields} | tr ' ' '\n' > {output.completion_flag}
        """

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/IBTrACS/0/wind_fields.txt
"""
