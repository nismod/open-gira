"""
Download STORM IBTrACS present climate synthetic tropical cyclone tracks and
tropical cyclone wind speed return periods

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085
"""


rule download_storm_basin_geometry:
    """
    Get a geometry file delineating the various storm basins
    """
    output:
        geometry = "{OUTPUT_DIR}/input/stormtracks/basins.kml"
    shell:
        """
        wget \
            https://data.4tu.nl/ndownloader/files/24060005 \
            --output-document {wildcards.OUTPUT_DIR}/input/stormtracks/basins.kml
        """


rule download_stormtracks_fixed:
    """
    Download storm return period and wind speed maps
    """
    output:
        "{OUTPUT_DIR}/input/stormtracks/fixed/STORM_FIXED_TC_WIND_SPEEDS_{REGION}.nc",
        "{OUTPUT_DIR}/input/stormtracks/fixed/STORM_FIXED_RETURN_PERIODS_{REGION}.nc"
    params:
        wind_speeds = "https://opendap.4tu.nl/thredds/fileServer/data2/uuid/779b9dfd-b0ff-4531-8833-aaa9c0cf6b5a/STORM_FIXED_TC_WIND_SPEEDS_{REGION}.nc",
        return_periods = "https://opendap.4tu.nl/thredds/fileServer/data2/uuid/779b9dfd-b0ff-4531-8833-aaa9c0cf6b5a/STORM_FIXED_RETURN_PERIODS_{REGION}.nc"
    shell:
        """
        wget \
            {params.wind_speeds} \
            {params.return_periods} \
            --directory-prefix={wildcards.OUTPUT_DIR}/input/stormtracks/fixed \
            --timestamping \
            --no-check-certificate
        """


rule download_stormtracks_events:
    """
    Download an archive of all storm event tracks for a given model (and some
    metadata, readmes, etc.)

    N.B. We rename the downloaded ZIP file from it's original name to
    archive.zip. This makes it easier to match on this file later. The mv
    command should fail if there's more than one zipfile as input (from the
    glob).
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/stormtracks/events/{STORM_MODEL}/archive.zip"
    shell:
        """
        wget \
            --input-file=config/hazard_resource_locations/storm_tracks_{wildcards.STORM_MODEL}.txt \
            --directory-prefix={wildcards.OUTPUT_DIR}/input/stormtracks/events/{wildcards.STORM_MODEL}/ \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        mv \
            {wildcards.OUTPUT_DIR}/input/stormtracks/events/{wildcards.STORM_MODEL}/*.zip \
            {wildcards.OUTPUT_DIR}/input/stormtracks/events/{wildcards.STORM_MODEL}/archive.zip
        """


rule extract_stormtracks_events:
    """
    Unzip a storm file for a basin we are interested in
    """
    input:
        rules.download_stormtracks_events.output.zip_file
    output:
        "{OUTPUT_DIR}/input/stormtracks/events/{STORM_MODEL}/{REGION}/{STORM_SAMPLE_BASENAME}.txt"
    shell:
        """
        unzip -o {input} {wildcards.STORM_SAMPLE_BASENAME}.txt \
            -d {wildcards.OUTPUT_DIR}/input/stormtracks/events/{wildcards.STORM_MODEL}/{wildcards.REGION}/
        """
