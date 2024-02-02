"""
Logic for handling exposure data.

Exposure data comprise the grid edge length `length_m`, exposed to wind speeds
in excess of a threshold.
"""

from open_gira.io import cached_json_file_read


rule aggregate_exposure_within_sample:
    """
    Take per-event exposure files with per-edge rows and aggregate into a per-edge file and a per-event file.
    """
    input:
        exposure_by_event = rules.electricity_grid_damages.output.exposure
    params:
        thresholds = config["transmission_windspeed_failure"]
    output:
        by_event = temp(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}_length_m_by_event.pq")),
        by_edge = temp(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}_length_m_by_edge.pq")),
    script:
        "./aggregate_grid_exposure.py"

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
        all_samples = protected("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/length_m_by_event.pq"),
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
        all_samples = protected("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/length_m_by_edge.pq"),
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
        "./plot_exposure_distributions.py"

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
        "./grid_exposure_by_admin_region.py"

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
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/EAE_{ADMIN_SLUG}.gpq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        ADMIN_SLUG=wildcards.ADMIN_SLUG  # str
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
        storm_set_exposure = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/EAE_{ADMIN_SLUG}.gpq"
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


def EAE_all_coarser_admin_levels(wildcards) -> list[str]:
    """
    Given an admin level, return a list of paths to the EAE file at that admin
    level and all the coarser levels down to and including level 0 (national).
    """
    max_level = int(wildcards.ADMIN_SLUG.split("-")[-1])
    return expand(
        "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/EAE_admin-level-{i}.gpq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        STORM_SET=wildcards.STORM_SET,
        i=range(max_level + 1),
    )

rule merge_exposure_admin_levels:
    """
    Take results at target admin level and gap fill with coarser admin levels
    where appropriate.
    """
    input:
        admin_levels = EAE_all_coarser_admin_levels,
    output:
        merged_admin_levels = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/EAE_{ADMIN_SLUG}-0.gpq"
    run:
        import logging

        import geopandas as gpd

        from open_gira.admin import merge_gadm_admin_levels

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        if len(input.admin_levels) == 1:
            # already at national level, no other region level data to merge
            gpd.read_parquet(input.admin_levels).to_parquet(output.merged_admin_levels)

        else:

            def read_and_label_ISO_A3(path: str) -> gpd.GeoDataFrame:
                """Requires a GID_n column with e.g. ZWE.9.3_1 or AFG.3_1 data."""
                df = gpd.read_parquet(path)
                GID_col, = df.columns[list(map(lambda s: s.startswith("GID_"), df.columns))]
                if "ISO_A3" in df.columns:
                    raise RuntimeError("Will not overwrite ISO_A3 column, quitting.")
                df["ISO_A3"] = df[GID_col].apply(lambda s: s.split(".")[0])
                return df

            # e.g. admin-level-3, admin-level-2, [admin-level-1, admin-level-0]
            # or admin-level-1, admin-level-0, []
            primary, secondary, *others = [read_and_label_ISO_A3(path) for path in sorted(input.admin_levels, reverse=True)]

            logging.info(f"Starting with: {set(primary.ISO_A3)}")
            merged = merge_gadm_admin_levels(primary, secondary)

            for other in others:
                merged = merge_gadm_admin_levels(merged, other)

            merged.reset_index(drop=True).sort_index(axis=1).to_parquet(output.merged_admin_levels)
            logging.info("Done")
