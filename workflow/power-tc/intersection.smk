"""
Estimate damages to electrical networks due to high wind speeds
"""

def storm_tracks_file_from_storm_set(wildcards) -> str:
    """
    Return a path to the processed tracks for a storm set.

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
    resources:
        mem_mb=60000
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
snakemake -c1 results/power/by_storm_set/IBTrACS/countries_impacted_by_storm.json
"""


rule electricity_grid_damages:
    """
    Degrade an electricity grid with a wind field, write out the exposure and
    disruption estimates.
    """
    input:
        # `threads_for_country` will fail unless this CSV is present when resources are set
        country_target_count=country_target_count_path,
        grid_splits = rules.rasterise_electricity_grid.output.geoparquet,
        wind_speeds = rules.estimate_wind_fields.output.wind_speeds,
        grid_edges = rules.create_power_network.output.edges,
        grid_nodes = rules.create_power_network.output.nodes,
    threads: threads_for_country
    resources:
        mem_mb = lambda wildcards: threads_for_country(wildcards) * 1_024 * 2.5
    output:
        exposure = protected(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/{STORM_SET}/{SAMPLE}/")),
        disruption = protected(directory("{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/{SAMPLE}/")),
    script:
        "./grid_disruption.py"

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/exposure/IBTrACS/0/
"""
