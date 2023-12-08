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
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}/{STORM_ID}.nc",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,  # str
        STORM_SET=wildcards.STORM_SET,  # str
        SAMPLE=wildcards.SAMPLE,  # str
        STORM_ID=storms  # list of str
    )

rule aggregate_exposure_within_sample:
    """
    Take per-event exposure files with per-edge rows and aggregate into a per-edge file and a per-event file.
    """
    input:
        exposure_by_event = exposure_by_storm_for_country_for_storm_set
    params:
        thresholds = config["transmission_windspeed_failure"]
    output:
        by_event = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}_length_m_by_event.pq"),
        by_edge = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}_length_m_by_edge.pq"),
    script:
        "../../../scripts/exposure/aggregate_grid_exposure.py"

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/0_length_m_by_event.pq
"""


def exposure_per_event_sample_files(wildcards) -> list[str]:
    """
    Return a list of paths, one for each sample.
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        rules.aggregate_exposure_within_sample.output.by_event,
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name]),
    )

rule aggregate_per_event_exposure_across_samples:
    """
    Take the per-sample exposure length files and combine them.
    """
    input:
        per_sample = exposure_per_event_sample_files,
    output:
        all_samples = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/length_m_by_event.pq",
    run:
        import pandas as pd

        df = pd.concat([pd.read_parquet(file_path) for file_path in input.per_sample])
        df.to_parquet(output.all_samples)

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/length_m_by_event.pq
"""


def exposure_per_edge_sample_files(wildcards) -> list[str]:
    """
    Return a list of paths, one for each sample.
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        rules.aggregate_exposure_within_sample.output.by_edge,
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=wildcards.COUNTRY_ISO_A3,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name]),
    )

rule aggregate_per_edge_exposure_across_samples:
    """
    Take the per-sample exposure length files and combine them.
    """
    input:
        per_sample = exposure_per_edge_sample_files,
    output:
        all_samples = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/length_m_by_edge.pq",
    run:
        import pandas as pd

        df = pd.concat([pd.read_parquet(file_path) for file_path in input.per_sample])
        df.groupby("edge").sum().to_parquet(output.all_samples)

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/length_m_by_edge.pq
"""


rule plot_event_exposure_distributions_for_country:
    """
    Plot the exposure distribution (over events) for a given country and storm set.
    """
    input:
        exposure_by_event = rules.aggregate_per_event_exposure_across_samples.output.all_samples
    output:
        country_event_distributions = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/length_m_event_dist/")
    script:
        "../../../scripts/exposure/plot_exposure_distributions.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/exposure/IBTrACS/length_m_event_dist
"""


rule exposure_by_admin_region:
    """
    Calculate expected annual exposure at given admin level.
    """
    input:
        tracks = storm_tracks_file_from_storm_set,  # storm dates (for time span -> expected annual exposure)
        exposure_by_event = rules.aggregate_per_event_exposure_across_samples.output.all_samples,  # event id list
        exposure_by_edge = rules.aggregate_per_edge_exposure_across_samples.output.all_samples,  # exposure data
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",  # regions to aggregate to
        grid_edges = rules.create_power_network.output.edges,  # edge linestrings
    output:
        expected_annual_exposure = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/EAE_{ADMIN_SLUG}.gpq",
    script:
        "../../../scripts/exposure/grid_exposure_by_admin_region.py"

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/exposure/IBTrACS/EAE_admin-level-1.gpq
"""


def exposure_summaries_for_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their summary files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/EAE_{ADMIN_LEVEL}.gpq",
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
        storm_set_exposure = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/EAE_{ADMIN_LEVEL}.gpq"
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
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/exposure/EAE_admin-level-2.gpq
"""