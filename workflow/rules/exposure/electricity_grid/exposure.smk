"""
Logic for handling exposure data.

Exposure data comprise the grid edge length `length_m`, exposed to wind speeds
in excess of a threshold.
"""

from open_gira.io import cached_json_file_read


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


rule concat_exposure_by_event:
    """
    Take per-event exposure files with per-edge rows and concatenate into a single file.
    """
    input:
        exposure_by_event = exposure_by_storm_for_country_for_storm_set,
    # read exposure files in parallel
    # use all available cores to block rest of execution (maximise memory available)
    threads: workflow.cores
    output:
        concatenated = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/exposure_by_event.parquet"
    script:
        "../../../scripts/exposure/concat_grid_exposure.py"

"""
Test with:
snakemake -n --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/exposure_by_event.parquet
"""


rule exposure_by_admin_region:
    """
    Calculate expected annual exposure at given admin level.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,
        exposure_by_edge_by_event = rules.concat_exposure_by_event.output.concatenated,
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        grid_edges = rules.create_power_network.output.edges,
    output:
        total_exposure_by_region = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{ADMIN_SLUG}.geoparquet",
    script:
        "../../../scripts/exposure/grid_exposure_by_admin_region.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/exposure/IBTrACS/admin-level-1.geoparquet
"""


rule plot_event_exposure_distributions_for_country:
    """
    Plot the exposure distribution (over events) for a given country and storm set.
    """
    input:
        exposure_by_event = rules.concat_exposure_by_event.output.concatenated
    output:
        country_event_distributions = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/event_distributions/")
    script:
        "../../../scripts/exposure/plot_exposure_distributions.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/exposure/IBTrACS/country_event_distribution
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

    Concatenates the regional summaries for expected annual exposure together.
    """
    input:
        exposure = exposure_summaries_for_storm_set
    output:
        storm_set_exposure = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/{ADMIN_LEVEL}.geoparquet"
    run:
        import geopandas as gpd
        import pandas as pd

        per_country_exposure = []
        for exposure_file in input.exposure:
            per_country_exposure.append(gpd.read_parquet(exposure_file))
        summary_file = gpd.GeoDataFrame(pd.concat(per_country_exposure)).reset_index(drop=True)
        summary_file.to_parquet(output.storm_set_exposure)

"""
Test with:
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/exposure/admin-level-2.geoparquet
"""
