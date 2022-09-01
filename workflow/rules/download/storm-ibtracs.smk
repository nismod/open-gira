"""
Download STORM IBTrACS present climate synthetic tropical cyclone tracks and
tropical cyclone wind speed return periods

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164
"""


rule download_stormtracks_fixed:
    """
    Download storm return period maps
    """
    output:
        STORMS_RETURN_PERIOD,
    shell:
        f"""
        wget \
            --input-file=config/hazard_resource_locations/storm_fixed_return.txt \
            --directory-prefix={config['output_dir']}/input/stormtracks/fixed \
            --timestamping \
            --no-check-certificate
        """


rule download_stormtracks_events:
    """
    Download an archive of all storm event tracks for a given model (and some
    metadata, readmes, etc.)
    """
    output:
        zip_file = os.path.join(f"{config['output_dir']}/input/stormtracks/events/", STORM_UNZIP_FILE)
    shell:
        f"""
        wget \
            --input-file=config/hazard_resource_locations/storm_tracks_{STORM_MODEL}.txt \
            --directory-prefix={config['output_dir']}/input/stormtracks/events \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        """


rule extract_stormtracks_events:
    """
    Unzip the storm files for basins we are interested in (config['regions'])
    """
    input:
        rules.download_stormtracks_events.output.zip_file
    output:
        STORMS_EVENTS
    run:
        # extract only the files in STORMS_EVENTS
        # STORMS_EVENTS is a list of full paths relative to repo root, use filenames instead
        concatenated_file_str = ' '.join(map(os.path.basename, STORMS_EVENTS))
        os.system(f"unzip -o {input} {concatenated_file_str} -d {config['output_dir']}/input/stormtracks/events")
