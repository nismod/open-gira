"""
Download STORM IBTrACS present climate synthetic tropical cyclone tracks and
tropical cyclone wind speed return period maps

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164/3
https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817/3
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085/4
https://data.4tu.nl/articles/dataset/STORM_Climate_Change_synthetic_tropical_cyclone_tracks/14237678
"""


rule download_storm:
    """
    Download an archive of all storm data for a given model (and some metadata,
    readmes, etc.). For event sets and return period maps.

    N.B. We rename the downloaded ZIP file from it's original name to
    archive.zip. This makes it easier to match on this file later. The mv
    command should fail if there's more than one zipfile as input (from the
    glob).
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/STORM/{EVENTS_OR_FIXED}/{STORM_MODEL}/archive.zip"
    shell:
        """
        if [[ "{wildcards.EVENTS_OR_FIXED}" == "events" ]]
        then
            RESOURCE_FILE="storm_tracks_{wildcards.STORM_MODEL}.txt"
        else
            RESOURCE_FILE="storm_fixed_{wildcards.STORM_MODEL}.txt"
        fi

        wget \
            --input-file=config/hazard_resource_locations/$RESOURCE_FILE \
            --directory-prefix={wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_FIXED}/{wildcards.STORM_MODEL}/ \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        mv \
            {wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_FIXED}/{wildcards.STORM_MODEL}/*.zip \
            {wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_FIXED}/{wildcards.STORM_MODEL}/archive.zip
        """

"""
Test with:
snakemake -c1 results/input/STORM/events/HadGEM-GC31-HM/archive.zip
"""

rule parse_storm:
    input:
        csv_dir="{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/raw"
    output:
        parquet="{OUTPUT_DIR}/storm_tracks/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet"
    script:
        "./parse_STORM.py"

"""
Test with:
snakemake -c1 results/storm_tracks/STORM-constant/0/tracks.geoparquet
"""


rule slice_storm:
    input:
        global_tracks="{OUTPUT_DIR}/storm_tracks/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "./slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/STORM-constant/0/tracks.geoparquet
"""

rule extract_storm_event:
    """
    Unzip a storm file for a basin we are interested in
    """
    input:
        "{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/archive.zip"
    output:
        "{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/{STORM_BASIN}/{STORM_SAMPLE_BASENAME}.txt"
    shell:
        """
        if [[ "{wildcards.STORM_MODEL}" == "constant" ]]
        then
            echo "Extracting current climate events"
            # version 4 has data nested in "VERSIE4" directory
            unzip -j -o {input} VERSIE4/{wildcards.STORM_SAMPLE_BASENAME}.txt \
                -d {wildcards.OUTPUT_DIR}/input/STORM/events/{wildcards.STORM_MODEL}/{wildcards.STORM_BASIN}/
        else
            echo "Extracting future climate events"
            unzip -o {input} {wildcards.STORM_SAMPLE_BASENAME}.txt \
                -d {wildcards.OUTPUT_DIR}/input/STORM/events/{wildcards.STORM_MODEL}/{wildcards.STORM_BASIN}/
        fi
        """

"""
Test with:
snakemake -c1 results/input/STORM/events/constant/NA/STORM_DATA_IBTRACS_NA_1000_YEARS_0.txt
snakemake -c1 results/input/STORM/events/HadGEM3-GC31-HM/NA/STORM_DATA_HadGEM3-GC31-HM_NA_1000_YEARS_0_IBTRACSDELTA.txt
"""


rule extract_all_storm_events:
    """
    Unzip all the storm files for a given model. Rename to appropriate file extension (CSV).
    """
    input:
        "{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/archive.zip"
    output:
        directory("{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/raw")
    shell:
        """
        MODEL_DIR={wildcards.OUTPUT_DIR}/input/STORM/events/{wildcards.STORM_MODEL}
        unzip -j -d $MODEL_DIR/raw {input}
        for FILE in $MODEL_DIR/raw/*.txt; do
            mv -- $FILE ${{FILE%.txt}}.csv
        done
        """

"""
Test with:
snakemake -c1 results/input/STORM/events/constant/raw/
"""


rule extract_storm_fixed_present:
    """
    Unzip a storm file
    """
    input:
        "{OUTPUT_DIR}/input/STORM/fixed/constant/archive.zip"
    output:
        "{OUTPUT_DIR}/input/STORM/fixed/constant/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    shell:
        """
        unzip -o {input} STORM_FIXED_RETURN_PERIODS_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif \
            -d {wildcards.OUTPUT_DIR}/input/STORM/fixed/constant/{wildcards.STORM_BASIN}/
        """

"""
Test with:
snakemake -c1 results/input/STORM/fixed/constant/NA/STORM_FIXED_RETURN_PERIODS_NA_20_YR_RP.tif
"""


rule rename_storm_fixed_present:
    input:
        rules.extract_storm_fixed_present.output
    output:
        "{OUTPUT_DIR}/input/STORM/fixed/constant/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_constant_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    shell:
        """
        mv {input} {output}
        """

"""
Test with:
snakemake -c1 results/input/STORM/fixed/constant/NA/STORM_FIXED_RETURN_PERIODS_constant_NA_20_YR_RP.tif
"""


rule extract_storm_fixed_future:
    """
    Unzip a storm file
    """
    input:
        "{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL_FUTURE}/archive.zip"
    output:
        "{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL_FUTURE}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL_FUTURE}_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    params:
        STORM_MODEL_UPPER = lambda wildcards: wildcards.STORM_MODEL_FUTURE.upper()
    shell:
        """
        # first try STORM_MODEL_FUTURE, otherwise fall back to trying the uppercase version
        unzip -o {input} STORM_FIXED_RETURN_PERIODS_{wildcards.STORM_MODEL_FUTURE}_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif \
            -d {wildcards.OUTPUT_DIR}/input/STORM/fixed/{wildcards.STORM_MODEL_FUTURE}/{wildcards.STORM_BASIN}/ \
        || \
        unzip -o {input} STORM_FIXED_RETURN_PERIODS_{params.STORM_MODEL_UPPER}_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif \
            -d {wildcards.OUTPUT_DIR}/input/STORM/fixed/{wildcards.STORM_MODEL_FUTURE}/{wildcards.STORM_BASIN}/ \
        """


rule wrap_storm_fixed:
    input:
        "{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    output:
        temp("{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_BASIN}_{STORM_RP}_YR_RP.wrapped.tif")
    shell:
        """
        gdalwarp -te -180 -60 180 60 {input} {output}
        """


rule mosaic_storm_fixed:
    """
    Merge basin return period maps to global extent
    """
    input:
        basin_tif_ep="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/EP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_EP_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_na="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/NA/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_NA_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_ni="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/NI/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_NI_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_si="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/SI/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_SI_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_sp="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/SP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_SP_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_wp="{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/WP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_WP_{STORM_RP}_YR_RP.wrapped.tif",
    output:
        "{OUTPUT_DIR}/input/STORM/fixed/{STORM_MODEL}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_RP}_YR_RP.tif"
    shell:
        """
        gdal_calc.py \
            -A {input.basin_tif_ep} \
            -B {input.basin_tif_na} \
            -C {input.basin_tif_ni} \
            -D {input.basin_tif_si} \
            -E {input.basin_tif_sp} \
            -F {input.basin_tif_wp} \
            --outfile={output} \
            --calc="numpy.max((A,B,C,D,E,F),axis=0)" \
            --NoDataValue=0
        """
