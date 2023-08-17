"""
Estimate damages to electrical networks due to high wind speeds
"""

from open_gira.io import cached_json_file_read


rule electricity_grid_damages:
    input:
        grid_splits = rules.rasterise_electricity_grid.output.geoparquet,
        wind_speeds = rules.estimate_wind_fields.output.wind_speeds,
        grid_edges = rules.create_power_network.output.edges,
        grid_nodes = rules.create_power_network.output.nodes,
    output:
        exposure = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc",
        disruption = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{STORM_ID}.nc"
    script:
        "../../scripts/intersect/grid_disruption.py"

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc
"""


def storm_tracks_file_from_storm_set(wildcards) -> str:
    """
    Return a path to the processed tracks for a stom set.

    Given e.g. IBTrACS_maria-2017, return results/input/storm-tracks/IBTrACS/tracks.geoparquet
    """
    storm_dataset = wildcards.STORM_SET.split("_")[0]
    return f"{wildcards.OUTPUT_DIR}/storm_tracks/{storm_dataset}/tracks.geoparquet"


checkpoint countries_intersecting_storm_set:
    """
    Find all countries which are likely to be affected by some storm set (as
    defined by a config["storm_sets"] JSON file).

    N.B. For a large set of tracks (> 1M) the spatial join can take a while.
    Currently the approach is to buffer the track points and not the countries.
    I have tried the inverse but it is much slower. Simplifying the country
    polygons does help (~10x faster), so that is implemented also.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,
        gadm_path = "{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
    output:
        country_set = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/countries_impacted.json",
        country_set_by_storm = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/countries_impacted_by_storm.json",
        storm_set_by_country = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storms_by_country_impacted.json",
    run:
        from collections import defaultdict
        import json
        import logging

        import geopandas as gpd
        import pandas as pd
        from tqdm import tqdm

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        # check for country intersections within some radius of track points
        search_radius_deg: float = config["max_track_search_radius_deg"]

        # read in track points
        logging.info("Reading tracks for storm set")
        tracks = gpd.read_parquet(input.tracks)

        # load storm ids of interest from file
        storm_set_path = config["storm_sets"][wildcards.STORM_SET]
        with open(storm_set_path, "r") as fp:
            storm_set = set(json.load(fp))

        if storm_set:
            logging.info(f"Filtering to {len(storm_set)} tracks of {wildcards.STORM_SET}")
            # subset points to the storm set of interest
            tracks_subset = tracks[tracks["track_id"].isin(storm_set)]
        else:
            # with no list of countries, assume we process all of them
            tracks_subset = tracks

        # only retain points where the observed max wind speed is greater than our threshold
        logging.info("Filter out points below minimum failure threshold")
        minimum_damage_threshold_ms = min(config["transmission_windspeed_failure"])
        tracks = tracks_subset.loc[tracks_subset.max_wind_speed_ms > minimum_damage_threshold_ms].copy()
        logging.info(f"Left with {len(tracks)} storm track points")
        logging.info(f"Buffer track points by {search_radius_deg} degrees for search")
        tracks.geometry = tracks.geometry.buffer(search_radius_deg)

        logging.info("Read country geometries")
        # combine natural earth (simple geom, 176 states) and GADM (detailed geom, 260 states)
        # for full country list, but without too much polygon detail for most countries
        # this speeds the later spatial join by ~10x
        gadm = gpd.read_parquet(input.gadm_path)[["GID_0", "geometry"]].rename(columns={"GID_0": "iso_a3"})
        natural_earth = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))[["iso_a3", "geometry"]]
        only_in_gadm_iso_a3 = list(set(gadm.iso_a3) - set(natural_earth.iso_a3))
        only_in_gadm = gadm.set_index("iso_a3", drop=True).loc[only_in_gadm_iso_a3].reset_index()
        only_in_gadm_copy = only_in_gadm.copy()
        # roughly simplify the remaining gadm geometries
        only_in_gadm_copy.geometry = only_in_gadm_copy.geometry.to_crs(epsg=3857).simplify(1000).to_crs(epsg=4326)
        countries = gpd.GeoDataFrame(pd.concat([natural_earth, only_in_gadm_copy]))

        # join track points to buffered country polygons
        logging.info("Perform spatial join of countries and filtered tracks")
        # current implementation takes ~4min to sjoin for 1,000 years of STORM-HadGEM
        intersecting_track_pts = pd.DataFrame(
            countries.sjoin(
                tracks[["max_wind_speed_ms", "geometry", "track_id"]],
                how="right"
            )[["track_id", "iso_a3"]]
        )

        logging.info("Parse join to identify which countries are impacted by which storms")
        country_set_by_storm = {}
        storm_set_by_country = defaultdict(list)
        for track_id, df in tqdm(intersecting_track_pts.groupby("track_id")):
            # unique ISO A3 country codes of likely affected countries
            if impacted_countries := list(set(df.loc[~df.iso_a3.isna(), "iso_a3"])):
                country_set_by_storm[track_id] = impacted_countries

                # build inverted map as we go
                for country in impacted_countries:
                    storm_set_by_country[country].append(track_id)

        logging.info("Write out to disk")

        # e.g. "AUS": ["200523501N234", "201201202S234"]
        with open(output.storm_set_by_country, "w") as fp:
            json.dump(dict(sorted(storm_set_by_country.items())), fp, indent=2)

        # e.g. "200523501N234": ["AUS", "PNG"]
        with open(output.country_set_by_storm, "w") as fp:
            json.dump(country_set_by_storm, fp, indent=2)

        # e.g. ["AUS", "PNG"]
        country_set = sorted(set().union(*country_set_by_storm.values()))
        with open(output.country_set, "w") as fp:
            json.dump(country_set, fp, indent=2)

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS_irma-2017/countries_impacted_by_storm.json
"""


def country_storm_paths_for_storm_set(wildcards):
    """
    Return list of paths of country-storm disruption
    """
    import json
    import os

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    country_set_by_storm = cached_json_file_read(json_file)

    paths = []
    for storm_id, countries in country_set_by_storm.items():
        paths.extend(
            expand(
                "results/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{STORM_ID}.nc",
                COUNTRY_ISO_A3=countries,  # list of str
                STORM_SET=wildcards.STORM_SET,  # str
                STORM_ID=storm_id  # str
            )
        )

    return paths


def country_storm_paths_for_storm(wildcards):
    """
    Given a STORM_ID and STORM_SET as a wildcard, lookup the countries that storm
    affects and return paths to their disruption files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    country_set_by_storm = cached_json_file_read(json_file)

    return expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{STORM_ID}.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=wildcards.STORM_ID  # str
    )


rule merge_countries_of_storm:
    """
    Merge disruption estimates from all countries a storm hit.
    """
    input:
        disruption = country_storm_paths_for_storm
    output:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/disruption_by_target.nc",
    run:
        import logging
        import os

        import xarray as xr

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Reading and pooling targets from all country datasets")
        pooled_targets = xr.concat([xr.open_dataset(path) for path in input.disruption], dim="target")

        # a few targets may have been processed under more than one country, keep the first instance
        logging.info("Dropping duplicates")
        pooled_targets = pooled_targets.drop_duplicates("target")

        # write to disk
        os.makedirs(os.path.dirname(output.by_target), exist_ok=True)
        logging.info("Writing pooled per-target disruption to disk")
        pooled_targets.to_netcdf(
            output.by_target,
            encoding={var: {"zlib": True, "complevel": 9} for var in pooled_targets.keys()}
        )

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017260N12310/disruption_by_target.nc
"""


def disruption_by_target_for_all_storms_in_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the storms in the set.

    Return a list of the disruption_by_target.nc file paths for every storm in the set.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    country_set_by_storm = cached_json_file_read(json_file)

    storms = list(country_set_by_storm.keys())

    return expand(
        "results/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/disruption_by_target.nc",
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=storms  # list of str
    )


rule disruption_by_storm_for_storm_set:
    """
    A target rule to generate the exposure and disruption netCDFs for all
    targets affected (across multiple countries) for each storm.
    """
    input:
        disruption = disruption_by_target_for_all_storms_in_storm_set
    output:
        completion_flag = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption.txt"
    shell:
        """
        # one output file per line
        echo {input.disruption} | tr ' ' '\n' > {output.completion_flag}
        """

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/IBTrACS/disruption.txt
"""


def exposure_by_storm_for_country_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the storms in the set impacting given COUNTRY_ISO_A3.

    Return a list of the relevant exposure netCDF file paths.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.storm_set_by_country
    storm_set_by_country = cached_json_file_read(json_file)

    storms = storm_set_by_country[wildcards.COUNTRY_ISO_A3]

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,  # str
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=storms  # list of str
    )


rule exposure_by_admin_region:
    """
    Calculate expected annual exposure at given admin level.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,
        exposure = exposure_by_storm_for_country_for_storm_set,
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        grid_edges = rules.create_power_network.output.edges,
    threads: 8  # read exposure files in parallel
    output:
        total_exposure_by_region = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{ADMIN_SLUG}.geoparquet",
        # TODO: per region event distributions
        # exposure_event_distribution_by_region = dir("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{ADMIN_SLUG}/")
    script:
        "../../scripts/exposure/grid_exposure_by_admin_region.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/exposure/IBTrACS/admin-level-1.geoparquet
"""


def exposure_summaries_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their summary files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{ADMIN_LEVEL}.geoparquet",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        ADMIN_LEVEL=wildcards.ADMIN_LEVEL  # str
    )


rule exposure_by_admin_region_for_storm_set:
    """
    A target rule to generate the exposure and disruption netCDFs for all
    targets affected (across multiple countries) for each storm.
    """
    input:
        disruption = exposure_summaries_for_storm_set
    output:
        storm_set_exposure = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/{ADMIN_LEVEL}.geoparquet"
    run:
        import geopandas as gpd
        import pandas as pd

        per_country_disruption = []
        for disruption_file in input.disruption:
            per_country_disruption.append(gpd.read_parquet(disruption_file))
        summary_file = gpd.GeoDataFrame(pd.concat(per_country_disruption)).reset_index(drop=True)
        summary_file.to_parquet(output.storm_set_exposure)

"""
Test with:
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/exposure/admin-level-2.geoparquet
"""
