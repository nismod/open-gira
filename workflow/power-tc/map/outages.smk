"""
Mapping electricity outages
"""


def storm_tracks(wildcards) -> str:
    """Return path to storm tracks"""
    # TODO: this function could also be used in storm_tracks.plot_storm_tracks

    # parent dataset of storm set e.g. IBTrACS for IBTrACS_maria-2017
    storm_dataset = wildcards.STORM_SET.split("_")[0]

    return f"{wildcards.OUTPUT_DIR}/storm_tracks/{storm_dataset}/tracks.geoparquet"


rule animate_electricity_outages:
    """
    Map electricity outages for all targets affected by a single storm.
    """
    input:
        targets = "{OUTPUT_DIR}/power/targets.geoparquet",
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/disruption_by_target.nc",
        tracks = storm_tracks,
    output:
        plot = directory("{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/outage_map")
    run:
        import geopandas as gpd
        import matplotlib
        from shapely.geometry import box
        import xarray as xr

        from open_gira.plot.outages import animate_outage_by_threshold

        # Use the Agg backend as we do not want to spawn a GUI:
        # UserWarning: Starting a Matplotlib GUI outside of the main thread will likely fail.
        matplotlib.use("Agg")

        # read in data
        global_targets = gpd.read_parquet(input.targets)
        disruption = xr.open_dataset(input.by_target)
        tracks = gpd.read_parquet(input.tracks)
        borders = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        event_id = wildcards.STORM_ID

        # subset tracks
        track = tracks[tracks.track_id == event_id]

        # select the affected targets
        targets = disruption.supply_factor.sel(
            event_id=event_id,
            threshold=min(disruption.threshold)
        ).to_pandas()
        non_null_targets = targets[~targets.isna()]
        undersupplied_targets = non_null_targets[non_null_targets < 1]

        # create a bounding box with a buffer
        aoi_targets = global_targets[global_targets.id.isin(undersupplied_targets.index.values)]
        aoi = box(*aoi_targets.geometry.centroid.total_bounds).buffer(3)

        animate_outage_by_threshold(event_id, output.plot, disruption.threshold.values, disruption, aoi, global_targets, borders, track)

"""
To test:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017260N12310/outage_map
"""
