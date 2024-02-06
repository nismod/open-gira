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
        "results/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}/{STORM_ID}.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        SAMPLE=wilcards.SAMPLE,  # str
        STORM_ID=wildcards.STORM_ID  # str
    )


rule disruption_merge_countries_of_storm:
    """
    Merge disruption estimates from all countries a storm hit.
    """
    input:
        disruption = country_storm_paths_for_storm
    output:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/{SAMPLE}/{STORM_ID}/disruption_by_target.nc",
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
        disruption_by_event = rules.electricity_grid_damages.output.disruption
    params:
        thresholds = config["transmission_windspeed_failure"]
    output:
        by_event = temp(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}_pop_affected_by_event.pq")),
        by_target = temp(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}_pop_affected_by_target.pq")),
    script:
        "./aggregate_grid_disruption.py"

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
        all_samples = protected("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_by_event.pq"),
    run:
        import pandas as pd

        df = pd.concat([pd.read_parquet(file_path) for file_path in input.per_sample])
        df.to_parquet(output.all_samples)

"""
Test with:
snakemake --cores 1 -- results/power/by_country/PRI/disruption/IBTrACS/pop_affected_by_event.pq
"""


rule disruption_pop_affected_return_periods:
    """
    Calculate how many people are affected by storms corresponding to given
    return periods. Also output the estimated recurrence time for the largest
    event of every year in the storm set.
    """
    input:
        per_event = rules.aggregate_per_event_disruption_across_samples.output,
        tracks = "{OUTPUT_DIR}/storm_tracks/{STORM_SET}/tracks.geoparquet",
    params:
        return_periods = config["return_period_years"]
    output:
        return_periods_raw = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_RP_raw.pq",
        return_periods_interpolated = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_RP.pq",
    run:
        import pandas as pd
        import numpy as np

        # read in the number of people affected by each event (event rows, threshold columns)
        pop_affected = pd.read_parquet(input.per_event)
        thresholds: np.ndarray = pop_affected.columns.values.copy()

        tracks = pd.read_parquet(input.tracks, columns=["track_id", "year"])
        tracks: pd.DataFrame = tracks.drop_duplicates("track_id").set_index("track_id")
        n_years: int = tracks.year.max() - tracks.year.min() + 1

        # bring in year data for the storms (common event_id index)
        pop_affected: pd.DataFrame = pop_affected.join(tracks)

        # year index, threshold columns, values are event_ids of max pop affected that year
        event_ids_of_annual_max: pd.DataFrame = pop_affected.groupby("year").idxmax()

        data_by_threshold: list[pd.DataFrame] = []
        for threshold in thresholds:

            # index into pop_affected data with event_ids of largest disruption per year, sort ascending
            df = pd.DataFrame(pop_affected.loc[event_ids_of_annual_max[threshold], threshold])

            # move speed threshold from column name into vales of new `threshold` column
            df = df.rename(columns={threshold: "pop_affected"})

            # add zero-valued entries for years with no record in the time span
            years_with_no_data = pd.DataFrame(
                index=pd.Index([f"no_data" for i in range(n_years - len(df))], name="event_id"),
                data={"pop_affected": np.zeros(n_years - len(df))}
            )
            # sort to ascending pop_affected
            df = pd.concat([df, years_with_no_data]).sort_values("pop_affected")

            # store threshold as data (will be set as index later)
            df["threshold"] = threshold

            # calculate recurrence intervals using 'plotting position' formula
            # most impactful event gets rank 1
            df["rank"] = range(n_years, 0, -1)

            # T = n + 1 / m, where T is expected value of return period, n is number of observations and m observation rank
            # for a derivation of this formula, see Gumbel 1958, §2.1.4 or Makkonen 2006
            # N.B. this is valid for any underlying distribution
            df["return_period_years"] = (n_years + 1) / df["rank"]

            # the Annual Exceedance Probability, AEP, is 1 / T

            df = df.reset_index().set_index(["threshold", "event_id"])

            data_by_threshold.append(df)

        # write out raw return periods as calculated for the largest event in each year
        data = pd.concat(data_by_threshold)
        data.to_parquet(output.return_periods_raw)

        # interpolate between the raw events to find population affected for a given set of return periods
        interpolated_RP_by_threshold = []
        for threshold in thresholds:
            df = pd.DataFrame(
                {
                    "threshold": threshold,
                    "return_period_years": params.return_periods,
                    # note that numpy interp expects the x grid, i.e.
                    # return_period_years data to be monotonically increasing
                    "pop_affected": np.interp(
                        params.return_periods,
                        data.xs(threshold, level="threshold")["return_period_years"],
                        data.xs(threshold, level="threshold")["pop_affected"]
                    )
                }
            )
            interpolated_RP_by_threshold.append(df)

        interpolated_RP = pd.concat(interpolated_RP_by_threshold).set_index(["threshold", "return_period_years"])
        interpolated_RP.to_parquet(output.return_periods_interpolated)

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/disruption/STORM-constant/pop_affected_RP.pq
"""


rule plot_pop_affected_return_periods:
    """
    Plot population affected as a function of return period.
    """
    input:
        return_periods = rules.disruption_pop_affected_return_periods.output.return_periods_raw
    params:
        best_estimate_threshold = config["best_estimate_windspeed_failure_threshold"]
    output:
        plot = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_RP.png",
    run:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        matplotlib.use("Agg")

        f, ax = plt.subplots(figsize=(12,6))

        data_by_threshold = pd.read_parquet(input.return_periods)
        thresholds = np.array(list(map(float, data_by_threshold.index.get_level_values(level="threshold").unique())))
        max_delta = np.max(np.abs(thresholds - params.best_estimate_threshold))

        cmap = matplotlib.colormaps["magma"]
        for threshold in thresholds:

            delta = np.abs(threshold - params.best_estimate_threshold)

            threshold_str = f"{threshold:.1f}"
            df = data_by_threshold.xs(threshold_str, level="threshold")

            if threshold == params.best_estimate_threshold:
                line_style = "--"
            else:
                line_style = "-"
            ax.step(
                df["return_period_years"],
                df["pop_affected"],
                label=f"{threshold_str}",
                color=cmap(delta / max_delta),
                ls=line_style,
            )

        ax.set_xscale("log")
        ax.set_xlabel("Return period [years]", labelpad=10)

        ax.set_ylabel(f"Population at risk of disconnection", labelpad=10)

        ax.grid(alpha=0.4, which="both")
        plt.subplots_adjust(right=0.8)
        ax.legend(
            title="Threshold [ms-1]",
            fontsize=8,
            bbox_to_anchor=(1.02, 1.0),
            loc='upper left'
        )
        ax.set_title(f"{wildcards.COUNTRY_ISO_A3}: {wildcards.STORM_SET}")

        f.savefig(output.plot)

"""
Test with:
snakemake -c1 -- results/power/by_country/PRI/disruption/STORM-constant/pop_affected_RP.png
"""


def pop_affected_return_periods_plot_path(wildcards) -> list[str]:
    """
    Return a list of paths, one for each country.
    """
    import json

    country_set_path = f"{wildcards.OUTPUT_DIR}/power/by_storm_set/{wildcards.STORM_SET}/countries_impacted.json"
    with open(country_set_path, "r") as fp:
        countries = json.load(fp)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_RP.png",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=countries,
        STORM_SET=wildcards.STORM_SET,
    )

rule plot_pop_affected_return_periods_storm_set:
    """"
    Target rule for population affected as a function of return period plots for
    all countries affected by STORM_SET.
    """
    input:
        plot_paths = pop_affected_return_periods_plot_path,
    output:
        flag = touch("{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/pop_affected_RP_plots.flag")

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/STORM-constant/pop_affected_RP_plots.flag
"""


def pop_affected_return_period_paths(wildcards) -> list[str]:
    """
    Return a list of paths, one for each country.
    """
    import json

    country_set_path = f"{wildcards.OUTPUT_DIR}/power/by_storm_set/{wildcards.STORM_SET}/countries_impacted.json"
    with open(country_set_path, "r") as fp:
        countries = json.load(fp)

    return expand(
        rules.disruption_pop_affected_return_periods.output.return_periods_interpolated,
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=countries,
        STORM_SET=wildcards.STORM_SET,
    )


rule disruption_pop_affected_return_period_map:
    """
    Number of people at risk of disconnection by country, by return period.
    """
    input:
        return_period_data = pop_affected_return_period_paths,
        admin_boundaries = "{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
    params:
        return_periods = config["return_period_years"]
    output:
        return_period_maps = directory("{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption/pop_affected_RP/"),
    run:
        import geopandas as gpd
        import pandas as pd

        countries = gpd.read_parquet(input.admin_boundaries)
        countries = countries.rename(columns={"GID_0": "ISO_A3"}).set_index("ISO_A3")

        data_by_country = []
        for path in input.return_period_data:
            iso_a3: str = path.split("/")[-4]
            df = pd.read_parquet(path)
            df["ISO_A3"] = iso_a3
            df = df.reset_index().set_index(["return_period_years", "threshold", "ISO_A3"])
            # threshold indicies are transformed into columns, long -> wide
            df = df.unstack("threshold")
            df.columns = df.columns.droplevel(0)
            data_by_country.append(df)

        data = pd.concat(data_by_country)

        os.makedirs(output.return_period_maps)
        for return_period in params.return_periods:
            # select one return period -- table now has country rows and threshold columns
            df = data.xs(return_period, level="return_period_years")

            # merge in geometry column, and name
            gdf = gpd.GeoDataFrame(df.join(countries))

            gdf.to_parquet(os.path.join(output.return_period_maps, f"{return_period}.gpq"))

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/STORM-constant/disruption/pop_affected_RP
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
        all_samples = protected("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_by_target.pq"),
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
        "./grid_disruption_by_admin_region.py"

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
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/EAPA_{ADMIN_SLUG}.gpq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,  # str
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET,  # str
        ADMIN_SLUG=wildcards.ADMIN_SLUG  # str
    )


rule disruption_by_admin_region_for_storm_set:
    """
    Concatenate the regional summaries for expected annual population affected.
    """
    input:
        disruption = disruption_summaries_for_storm_set
    output:
        storm_set_disruption = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption/EAPA_{ADMIN_SLUG}.gpq"
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
