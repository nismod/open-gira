"""
Estimate max wind speed at infrastructure asset locations per event
"""

checkpoint intersect_winds_indiv:
    """
    Find windspeeds for every storm for each grid cell.
    """
    conda: "../../../environment.yml"
    input:
        metadata = rules.world_splitter.output.global_metadata,
        edges_split = "{OUTPUT_DIR}/processed/power/{BOX}/edges_split_{BOX}.geoparquet",
        storm_files = STORM_EVENTS,  # TODO limited to specific storms if any
        # Do we want ws output files per group of storms? (basin/model/sample)
        # TODO deal with inconsistency in present vs future filename patterns
        # storm_file="{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/{STORM_BASIN}/STORM_DATA_{STORM_GCM_OR_IBTRACS}_{STORM_BASIN}_1000_YEARS_{SAMPLE_IBTRACSDELTA}.txt",
    output:
        # TODO figure out if this should be per storm or per sample?
        # TODO BASIN and BOX is redundant - can we skip BASIN? would do some lookup function thing
        # STORM_SUBSET - key in config dict
        wind_speeds = "{OUTPUT_DIR}/processed/power/exposure/{STORM_SUBSET}/{STORM_MODEL}/{BASIN}/ws_{BOX}_{SAMPLE}.parquet",
        # TODO gather into "{OUTPUT_DIR}/processed/power/exposure/{STORM_SUBSET}.parquet",
    script:
        "../../scripts/intersect/intersect_3_wind_extracter.py"
