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


checkpoint electricity_grid_damages:
    input:
        grid_splits = rules.rasterise_electricity_grid.output.geoparquet,
        wind_speeds = rules.estimate_wind_fields.output.wind_speeds,
        grid_edges = rules.create_power_network.output.edges,
        grid_nodes = rules.create_power_network.output.nodes,
    threads:
        config["processes_per_parallel_job"]
    output:
        damages = directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}"),
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


def countries_endangered_by_storm_set(wildcards):
    """
    Return list of paths of country exposure directories to generate
    """
    import json
    import os

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    with open(json_file, "r") as fp:
        country_set = json.load(fp)

    file_paths = expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}",
        COUNTRY_ISO_A3=country_set,
        STORM_SET=wildcards.STORM_SET
    )

    return file_paths


rule storm_set_damages:
    """
    A target rule for storm sets, expanding to calculate exposure for all
    affected countries.

    To indicate completion, copy storm set JSON file to results directory.
    """
    input:
        countries_endangered_by_storm_set
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


# some example nested globs/checkpoints
#
#   def aggregate(wildcards):
#       outputs_i = glob.glob(f"{checkpoints.first.get().output}/*/")
#       outputs_i = [output.split('/')[-2] for output in outputs_i]
#       split_files = []
#       for i in outputs_i:
#           outputs_j = glob.glob(f"{checkpoints.second.get(i=i).output}/*/")
#           outputs_j = [output.split('/')[-2] for output in outputs_j]
#           for j in outputs_j:
#               split_files.append(f"copy/{i}/{j}/test2.txt")
#
#       return split_files
#
#   def aggregate_decompress_plass(wildcards):
#       checkpoint_output = checkpoints.decompress_plass.get(**wildcards).output[0]
#       file_names = expand(
#           "outputs/cd-hit95/{mag}.cdhit95.faa.bwt",
#           mag = glob_wildcards(
#               os.path.join(
#                   checkpoint_output,
#                   "{mag}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup"
#               )
#           ).mag
#       )
#       return file_names


# this is a lead, sorta
# works if you have a checkpoint get call before the (actually useful) glob
# this delays execution of this function by snakemake until electricity_grid_damages has completed
# alas electricity_grid_damages needs a COUNTRY_ISO_A3 wildcard to run
# here we give it a constant (but relevant) one to make it run
#
#   def countries_hit_by_storm(wildcards):
#       """
#       Return list of paths of per-country exposure files for a given storm.
#       """
#       # hack: use checkpoint.get to trigger delayed execution of this function
#       checkpoints.electricity_grid_damages.get(**wildcards, COUNTRY_ISO_A3="PRI").output.damages
#
#       from glob import glob
#       exposure_paths = glob(f"{wildcards.OUTPUT_DIR}/power/by_country/*/exposure/{wildcards.STORM_SET}/{wildcards.STORM_ID}.nc")
#
#       return exposure_paths

# modification of existing input function for storm sets, works in that it generates a list of netCDFs to produce
# broken in that it doesn't trigger generation of the upstream (electricity_grid_damages does not run as a result)
# list of countries also the superset of what actually gets produced
#
def countries_endangered_by_storm(wildcards):
    """
    Return list of paths of country exposure directories to generate
    """
    import json

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set_by_storm
    with open(json_file, "r") as fp:
        country_set_by_storm = json.load(fp)

    file_paths = expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{STORM_ID}.nc",
        COUNTRY_ISO_A3=country_set_by_storm[wildcards.STORM_ID],  # list of ISO strs
        STORM_SET=wildcards.STORM_SET,
        STORM_ID=wildcards.STORM_ID
    )

    # some of the countries we thought might be effected are actually unscathed
    # so output files aren't necessarily created for all the countries we expected
    # filter out these missing file paths
    # this is probably dubious snakemake usage, so if there's a better solution, please change!
    return list(filter(os.path.exists, file_paths))


rule concat_countries_of_storm:
    """
    Concatenate exposure estimates from all countries a storm hit.
    """
    input:
        # require grid simulation to have completed first
        completion_flag = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/storm_set.json",
        per_country_exposure = countries_endangered_by_storm
    output:
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/{STORM_ID}.nc",
    run:
        import logging

        import xarray as xr

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Reading and pooling targets from all country datasets")
        pooled_targets = xr.concat([xr.open_dataset(path) for path in input.per_country_exposure], dim="target")

        # a few targets may have been processed under more than one country, keep the first instance
        logging.info("Dropping duplicates")
        pooled_targets = pooled_targets.drop_duplicates("target")

        # write to disk
        logging.info("Writing pooled per-target exposure to disk")
        pooled_targets.to_netcdf(
            output.by_target,
            encoding={var: {"zlib": True, "complevel": 9} for var in pooled_targets.keys()}
        )

"""
Test with:
snakemake -c1 results/power/by_storm_set/IBTrACS_maria-2017/2017260N12310.nc"
"""


def concat_storm_filepaths(wildcards) -> list[str]:
    """
    Return a list of exposure storm file paths.
    """
    storm_set_path = lambda wildcards: config["storm_sets"][wildcards.STORM_SET]
    with open(storm_set_path, "r") as fp:
        storm_ids_of_storm_set = json.load(fp)
    ret = expand(
        "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/{STORM_ID}.nc",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        STORM_SET=wildcards.STORM_SET,
        STORM_ID=storm_ids_of_storm_set
    )
    print(ret)
    return ret


rule concat_countries_of_storm_set:
    """
    Concatenate across countries to individual storm files for all storms in a storm set.
    """

    input:
        storm_files = concat_storm_filepaths
    output:
        "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/done.txt"
    shell:
        """
        touch {output}
        """
