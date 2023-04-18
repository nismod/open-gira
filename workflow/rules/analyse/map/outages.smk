"""
Mapping electricity outages
"""


def storm_tracks(wildcards) -> str:
    """Return path to storm tracks"""
    # TODO: this function could also be used in storm_tracks.plot_storm_tracks

    # parent dataset of storm set e.g. IBTrACS for IBTrACS_maria-2017
    storm_dataset = wildcards.STORM_SET.split("_")[0]

    # TODO: these folders should really have a more consistent structure
    # e.g., move all the processed stuff to results/power/events/{storm_dataset}/tracks.geoparquet
    # then we can remove this conditional
    if storm_dataset == "IBTrACS":
        return f"{wildcards.OUTPUT_DIR}/input/IBTrACS/processed/v4.geoparquet"
    elif storm_dataset == "STORM-constant":
        return f"{wildcards.OUTPUT_DIR}/input/STORM/events/STORM-constant/processed.geoparquet"
    else:
        raise NotImplementedError(f"Do not recognise {storm_dataset=}, please add to rule if desired.")


rule animate_electricity_outages:
    """
    Map electricity outages for all targets affected by a single storm.
    """
    input:
        targets = "{OUTPUT_DIR}/power/targets.geoparquet",
        by_target = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/by_storm/{STORM_ID}/exposure_by_target.nc",
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
        targets = gpd.read_parquet(input.targets)
        exposure = xr.open_dataset(input.by_target)
        tracks = gpd.read_parquet(input.tracks)
        borders = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        event_id = wildcards.STORM_ID

        # subset tracks
        track = tracks[tracks.track_id == event_id]

        # find all targets within X degrees of track
        track_bbox = box(*track.geometry.total_bounds)
        aoi = track_bbox.buffer(5)
        aoi_targets = targets[targets.within(aoi)]

        # shrink AOI to ignore track in the open ocean
        track_target_aoi = box(*aoi_targets.geometry.centroid.total_bounds)

        animate_outage_by_threshold(event_id, output.plot, exposure.threshold.values, exposure, track_target_aoi, aoi_targets, borders, track)

"""
To test:
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017260N12310/outage_map
"""
