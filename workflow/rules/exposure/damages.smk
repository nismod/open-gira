"""
Estimate damage fraction for exposed assets
Damage fraction is a function of hazard intensity (expressed as damage curves)
"""


rule direct_damages:
    input:
        unsplit = rules.create_transport_network.output.edges,  # for pre-intersection geometry
        exposure = rules.rasterise_osm_network.output.geoparquet
    output:
        damage_fraction = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/fraction_per_RP/{SLICE_SLUG}.geoparquet",
        damage_cost = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/cost_per_RP/{SLICE_SLUG}.geoparquet",
        expected_annual_damages = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD/{SLICE_SLUG}.geoparquet",
        return_period_and_ead = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_RP/{SLICE_SLUG}.geoparquet",
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        # determine the hazard type from the hazard slug, e.g. flood, earthquake, storm
        hazard_type=lambda wildcards: config["hazard_types"][wildcards.HAZARD_SLUG.replace('hazard-', '')]
    script:
        "../../scripts/direct_damages.py"

"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/slice-5.geoparquet
"""


rule plot_damage_distributions:
    input:
        damages = "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction_per_RP.geoparquet"
    output:
        plots = directory("{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction_plots")
    script:
        "../../scripts/plot_damage_distributions.py"

"""
Test with:
snakemake --cores 1 results/egypt-latest_filter-road/hazard-aqueduct-river/damage_fraction_plots
"""


rule electricity_grid_damages:
    input:
        grid_splits = rules.rasterise_electricity_grid.output.geoparquet,
        wind_speeds = rules.estimate_wind_fields.output.wind_speeds,
        grid_edges = rules.create_power_network.output.edges,
        grid_nodes = rules.create_power_network.output.nodes,
    threads:
        config["processes_per_parallel_job"]
    output:
        damages = "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}.nc",
    script:
        "../../scripts/intersect/grid_disruption.py"

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/exposure/IBTrACS.nc
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
        country_set = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/countries_hit.json"
    run:
        import json

        import pandas as pd
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

        danger_zone = ibtracs_subset.copy()
        danger_zone.geometry = ibtracs_subset.geometry.buffer(point_buffer_deg)

        countries = gpd.read_parquet(input.admin_bounds).rename(columns={"GID_0": "iso_a3"})

        # join buffered points to countries
        hit = countries.sjoin(danger_zone[["max_wind_speed_ms", "geometry"]], how="right")

        # only retain points where the observed max wind speed is greater than our threshold
        minimum_damage_threshold_ms = min(config["transmission_windspeed_failure"])
        hit = hit.loc[hit.max_wind_speed_ms > minimum_damage_threshold_ms, "iso_a3"]

        # unique ISO A3 country codes of likely affected countries
        hit_iso_a3: list[str] = list(set(hit[~hit.isna()].values))

        with open(output.country_set, "w") as fp:
            json.dump(sorted(hit_iso_a3), fp, indent=2)

"""
Test with:
snakemake -c1 results/power/by_storm_set/black_marble_validation/countries_hit.json
"""


def at_risk_countries(wildcards):
    """
    Return list of paths of country exposure files to generate
    """
    import json

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    with open(json_file, "r") as fp:
        country_set = json.load(fp)

    return expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}.nc",
        COUNTRY_ISO_A3=country_set,
        STORM_SET=wildcards.STORM_SET
    )


rule storm_set_damages:
    """
    A target rule for storm sets, expanding to calculate exposure for all
    affected countries.

    To indicate completion, copy storm set JSON file to results directory.
    """
    input:
        at_risk_countries
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
snakemake -c1 results/power/by_storm_set/black_marble_validation/storm_set.json
"""


rule combine_storm_set_exposure:
    """
    Concatenate per-country exposure results and save to disk. Additionally,
    aggregate to country level and save to disk.
    """
    input:
        per_country_exposure = at_risk_countries
    output:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure_by_target.nc",
        by_country = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure_by_country.nc"
    run:
        from collections import defaultdict

        import pandas as pd
        import xarray as xr


        # ISO -> country data
        datasets: dict[str: xr.Dataset] = {
            path.split("/")[3]:
            xr.open_dataset(path) for path in input.per_country_exposure
        }

        # update target ids to make them globally unique
        # e.g. USA-1173 or PRI-101
        unique_target_datasets = []
        for iso, ds in datasets.items():
            # our empty datasets don't have coordinates
            if hasattr(ds, "target") and hasattr(ds, "event_id"):
                unique_target_datasets.append(
                    ds.assign_coords(
                        {
                            "target":
                            [f"{iso}-{target_id}" for target_id in ds.target.values]
                        }
                    )
                )

        # combine country datasets of targets to pool all targets together in one file
        # N.B. this is a quite sparse, with every target from every country for every storm in the set
        # for IBTrACS, there are ~800M values per variable, and both variables are 96% NaN
        concat = xr.concat(unique_target_datasets, dim="target").sortby("event_id")
        # write to disk
        concat.to_netcdf(output.by_target, encoding={var: {"zlib": True, "complevel": 9} for var in concat.keys()})

        # mapping from country to list of target ids in country
        targets_by_country: dict[str, list[str]] = defaultdict(list)
        for uid in concat.target.values:
            targets_by_country[uid.split("-")[0]].append(uid)

        # sum data across targets to country level
        countries = list(datasets.keys())
        country_data = []
        for country in countries:
            country_data.append(
                xr.DataArray(
                    concat.customers_affected.sel(dict(target=targets_by_country[country])).sum(dim="target")
                )
            )
        dims = ["country", "event_id", "threshold"]
        by_country: xr.DataArray = xr.concat(country_data, dim=pd.Index(list(countries), name="country"))
        by_country.to_netcdf(output.by_country, encoding={"customers_affected": {"zlib": True, "complevel": 9}})

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS_maria-2017/exposure_by_country.nc"
"""
