"""
Logic for handling disruption data.

Disruption data comprise estimates of `supply_factor`, the change in available
power supply to an electricity consuming target.
"""

from open_gira.io import cached_json_file_read


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
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/0/2017260N12310/disruption_by_target.nc
"""


rule aggregate_disruption_within_sample:
    """
    Take per-event disruption files with per-target rows (for all of a storm set
    sample) and aggregate into a per-target file and a per-event file.
    """
    input:
        disruption_by_event = disruption_by_storm_for_country_for_storm_set
    params:
        thresholds = config["transmission_windspeed_failure"]
    output:
        by_event = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}_pop_affected_by_event.pq"),
        by_target = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}_pop_affected_by_target.pq"),
    script:
        "../../../scripts/exposure/aggregate_grid_disruption.py"

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/disruption/IBTrACS/0_pop_affected_by_event.pq
"""


def disruption_per_event_sample_files(wildcards) -> list[str]:
    """
    Return a list of paths, one for each sample.
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        rules.aggregate_disruption_within_sample.output.by_event,
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name]),
    )

rule aggregate_per_event_disruption_across_samples:
    """
    Take the per-sample customers affected files and combine them.
    """
    input:
        per_sample = disruption_per_event_sample_files,
    output:
        all_samples = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_by_event.pq",
    run:
        import pandas as pd

        df = pd.concat([pd.read_parquet(file_path) for file_path in input.per_sample])
        df.to_parquet(output.all_samples)

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/disruption/IBTrACS/pop_affected_by_event.pq
"""


def disruption_per_target_sample_files(wildcards) -> list[str]:
    """
    Return a list of paths, one for each sample.
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        rules.aggregate_disruption_within_sample.output.by_target,
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name]),
    )

rule aggregate_per_target_disruption_across_samples:
    """
    Take the per-sample customers affected files and combine them.
    """
    input:
        per_sample = disruption_per_target_sample_files,
    output:
        all_samples = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_by_target.pq",
    run:
        import pandas as pd

        df = pd.concat([pd.read_parquet(file_path) for file_path in input.per_sample])
        df.groupby("target").sum().to_parquet(output.all_samples)

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/disruption/IBTrACS/pop_affected_by_target.pq
"""


rule disruption_by_admin_region:
    """
    Calculate expected annual population affected at given admin level.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,
        disruption_by_target = rules.aggregate_per_target_disruption_across_samples.output.all_samples,
        disruption_by_event = rules.aggregate_per_event_disruption_across_samples.output.all_samples,
        targets = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
    output:
        expected_annual_disruption = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/EAPA_{ADMIN_SLUG}.gpq",
    script:
        "../../../scripts/exposure/grid_disruption_by_admin_region.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/disruption/IBTrACS/EAPA_admin-level-1.gpq
"""


def disruption_summaries_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their summary files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/EAPA_{ADMIN_LEVEL}.gpq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        ADMIN_LEVEL=wildcards.ADMIN_LEVEL  # str
    )


rule disruption_by_admin_region_for_storm_set:
    """
    Concatenate the regional summaries for expected annual population affected.
    """
    input:
        disruption = disruption_summaries_for_storm_set
    output:
        storm_set_disruption = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption/EAPA_{ADMIN_LEVEL}.gpq"
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
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/disruption/EAPA_admin-level-2.gpq
"""
