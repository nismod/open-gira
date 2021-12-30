"""Download STORM IBTrACS present climate synthetic tropical cyclone tracks and tropical cyclone 
wind speed return periods

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164
"""

CYCLONE_REGIONS = ["EP", "NA", "NI", "SI", "SP", "WP"]

out_fixed = expand(
    os.path.join(DATA_DIR, "stormtracks", "fixed", "STORM_FIXED_{param}_{region}.nc"),
    region=CYCLONE_REGIONS,
    param=["RETURN_PERIODS", "TC_WIND_SPEEDS"],
)

out_events = expand(
    os.path.join(
        DATA_DIR,
        "stormtracks",
        "events",
        "STORM_DATA_IBTRACS_{region}_1000_YEARS_{num}.txt",
    ),
    region=CYCLONE_REGIONS,
    num=list(range(0, 10)),
)


rule download_stormtracks_fixed:
    output:
        out_fixed,
    shell:
        """
        mkdir -p data/stormtracks/fixed/extracted
        wget \
            --input-file=workflow/scripts/storm_fixed_return.txt \
            --directory-prefix=data/stormtracks/fixed \
            --timestamping \
            --no-check-certificate
        """


rule download_stormtracks_events:
    output:
        out_events,
    shell:
        """
        wget \
            --input-file=workflow/scripts/storm_tracks.txt \
            --directory-prefix=data/stormtracks/events \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        unzip -o data/stormtracks/events/STORM_DATA3.zip -d data/stormtracks/events
        """
