"""
Estimate damages to electrical networks due to high wind speeds
"""


rule electricity_grid_damages:
    input:
        grid_splits = rules.rasterise_electricity_grid.output.geoparquet,
        wind_speeds = rules.estimate_wind_fields.output.wind_speeds,
        grid_edges = rules.create_power_network.output.edges,
        grid_nodes = rules.create_power_network.output.nodes,
    output:
        exposure = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc"
    script:
        "../../scripts/intersect/grid_disruption.py"

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc
"""


checkpoint countries_intersecting_storm_set:
    """
    Find all countries which are likely to be affected by some storm set (as
    defined by a config["storm_sets"] JSON file).
    """
    input:
        ibtracs = "{OUTPUT_DIR}/input/IBTrACS/processed/v4.geoparquet",
        admin_bounds = "{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
    output:
        country_set = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/countries_impacted.json",
        country_set_by_storm = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/countries_impacted_by_storm.json",
        storm_set_by_country = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storms_by_country_impacted.json",
    run:
        from collections import defaultdict
        import json

        import geopandas as gpd

        # check for country intersections within some radius of track points
        point_buffer_deg: float = config["max_track_search_radius_deg"]

        # read in track points
        ibtracs = gpd.read_parquet(input.ibtracs)

        # load storm ids of interest from file
        storm_set_path = config["storm_sets"][wildcards.STORM_SET]
        with open(storm_set_path, "r") as fp:
            storm_set = set(json.load(fp))

        if storm_set:
            # subset points to the storm set of interest
            ibtracs_subset = ibtracs[ibtracs["track_id"].isin(storm_set)]
        else:
            # with no list of countries, assume we process all of them
            ibtracs_subset = ibtracs

        tracks = ibtracs_subset.copy()
        tracks.geometry = tracks.geometry.buffer(point_buffer_deg)

        countries = gpd.read_parquet(input.admin_bounds).rename(columns={"GID_0": "iso_a3"})

        # join buffered track points to country IDs
        intersection = countries.sjoin(
            tracks[["max_wind_speed_ms", "geometry", "track_id"]],
            how="right"
        )

        # only retain points where the observed max wind speed is greater than our threshold
        minimum_damage_threshold_ms = min(config["transmission_windspeed_failure"])
        above_threshold_intersection = intersection.loc[
            intersection.max_wind_speed_ms > minimum_damage_threshold_ms,
            ["iso_a3", "track_id"]
        ]

        country_set_by_storm = {}
        storm_set_by_country = defaultdict(list)
        for track_id, df in above_threshold_intersection.groupby("track_id"):
            # unique ISO A3 country codes of likely affected countries
            if impacted_countries := list(set(df.loc[~df.iso_a3.isna(), "iso_a3"])):
                country_set_by_storm[track_id] = impacted_countries

                # build inverted map as we go
                for country in impacted_countries:
                    storm_set_by_country[country].append(track_id)

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
    Return list of paths of country-storm exposure
    """
    import json
    import os

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    with open(json_file, "r") as fp:
        country_set_by_storm = json.load(fp)

    paths = []
    for storm_id, countries in country_set_by_storm.items():
        paths.extend(
            expand(
                "results/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc",
                COUNTRY_ISO_A3=countries,  # list of str
                STORM_SET=wildcards.STORM_SET,  # str
                STORM_ID=storm_id  # str
            )
        )

    return paths


rule storm_set_damages:
    """
    A target rule for storm sets, expanding to calculate exposure for all
    affected countries.

    To indicate completion, copy storm set JSON file to results directory.
    """
    input:
        country_storm_paths_for_storm_set
    params:
        storm_set_path = lambda wildcards: config["storm_sets"][wildcards.STORM_SET]
    output:
        completion_flag = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_set.json"
    shell:
        """
        cp {params.storm_set_path} {output.completion_flag}
        """

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS/storm_set.json
"""


def country_storm_paths_for_storm(wildcards):
    """
    Given a STORM_ID and STORM_SET as a wildcard, lookup the countries that storm
    affects and return paths to their exposure files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    with open(json_file, "r") as fp:
        country_set_by_storm = json.load(fp)

    return expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=wildcards.STORM_ID  # str
    )


rule merge_countries_of_storm:
    """
    Merge exposure estimates from all countries a storm hit.
    """
    input:
        exposure = country_storm_paths_for_storm
    output:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/exposure_by_target.nc",
    run:
        import logging
        import os

        import xarray as xr

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Reading and pooling targets from all country datasets")
        pooled_targets = xr.concat([xr.open_dataset(path) for path in input.exposure], dim="target")

        # a few targets may have been processed under more than one country, keep the first instance
        logging.info("Dropping duplicates")
        pooled_targets = pooled_targets.drop_duplicates("target")

        # write to disk
        os.makedirs(os.path.dirname(output.by_target), exist_ok=True)
        logging.info("Writing pooled per-target exposure to disk")
        pooled_targets.to_netcdf(
            output.by_target,
            encoding={var: {"zlib": True, "complevel": 9} for var in pooled_targets.keys()}
        )

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017260N12310/exposure_by_target.nc"
"""
