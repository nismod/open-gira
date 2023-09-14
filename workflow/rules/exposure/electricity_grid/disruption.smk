"""
Logic for handling disruption data.

Disruption data comprise estimates of `supply_factor`, the change in available
power supply to an electricity consuming target.
"""

from open_gira.io import cached_json_file_read


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
                "results/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/by_storm/{STORM_ID}.nc",
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
        "results/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/by_storm/{STORM_ID}.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=wildcards.STORM_ID  # str
    )


rule disruption_merge_countries_of_storm:
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
        pooled_targets.to_netcdf(output.by_target)

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


def disruption_by_storm_for_country_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the storms in the set impacting given COUNTRY_ISO_A3.

    Return a list of the relevant disruption netCDF file paths.
    """
    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.storm_set_by_country
    storm_set_by_country = cached_json_file_read(json_file)

    storms = storm_set_by_country[wildcards.COUNTRY_ISO_A3]

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/by_storm/{STORM_ID}.nc",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,  # str
        STORM_SET=wildcards.STORM_SET,  # str
        STORM_ID=storms  # list of str
    )


rule disruption_by_admin_region:
    """
    Calculate expected annual population affected at given admin level.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,
        disruption = disruption_by_storm_for_country_for_storm_set,
        targets = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
    threads: 8  # read exposure files in parallel
    output:
        total_disruption_by_region = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{ADMIN_SLUG}.geoparquet",
        # TODO: per region event distributions
        # disruption_event_distribution_by_region = dir("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{ADMIN_SLUG}/")
    script:
        "../../../scripts/exposure/grid_disruption_by_admin_region.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/disruption/IBTrACS/admin-level-1.geoparquet
"""


def disruption_summaries_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their summary files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{ADMIN_LEVEL}.geoparquet",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        ADMIN_LEVEL=wildcards.ADMIN_LEVEL  # str
    )


rule disruption_by_admin_region_for_storm_set:
    """
    A target rule to generate the exposure and disruption netCDFs for all
    targets affected (across multiple countries) for each storm.

    Concatenates the regional summaries for expected annual population disruption together.
    """
    input:
        disruption = disruption_summaries_for_storm_set
    output:
        storm_set_disruption = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption/{ADMIN_LEVEL}.geoparquet"
    run:
        import geopandas as gpd
        import pandas as pd

        per_country_disruption = []
        for disruption_file in input.disruption:
            per_country_disruption.append(gpd.read_parquet(disruption_file))
        summary_file = gpd.GeoDataFrame(pd.concat(per_country_disruption)).reset_index(drop=True)
        summary_file.to_parquet(output.storm_set_disruption)

"""
Test with:
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/disruption/admin-level-2.geoparquet
"""
