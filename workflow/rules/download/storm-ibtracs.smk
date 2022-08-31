"""Download STORM IBTrACS present climate synthetic tropical cyclone tracks and tropical cyclone 
wind speed return periods

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164
"""


rule download_stormtracks_fixed:
    output:
        STORMS_RETURN_PERIOD,
    shell:
        f"""
        wget \
            --input-file=workflow/scripts/storm_fixed_return.txt \
            --directory-prefix={config['output_dir']}/input/stormtracks/fixed \
            --timestamping \
            --no-check-certificate
        """


rule download_stormtracks_events:
    output:
        STORMS_EVENTS,
    shell:
        f"""
        wget \
            --input-file=workflow/scripts/storm_tracks_{STORM_MODEL}.txt \
            --directory-prefix={config['output_dir']}/input/stormtracks/events \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        unzip -o {config['output_dir']}/input/stormtracks/events/{UNZIP_FILE} -d {config['output_dir']}/input/stormtracks/events
        """
