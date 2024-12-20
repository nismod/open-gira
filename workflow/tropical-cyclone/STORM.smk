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


rule download_STORM:
    """
    Download an archive STORM wind speed data for a given model (and some
    metadata, readmes, etc.). For event sets (tracks) and raster maps (fixed
    RP, or fixed wind speed).

    N.B. We rename the downloaded ZIP file from it's original name to
    archive.zip. This makes it easier to match on this file later. The mv
    command should fail if there's more than one zipfile as input (from the
    glob).

    Test with:
    snakemake -c1 results/input/STORM/events/HadGEM3-GC31-HM/archive.zip
    snakemake -c1 results/input/STORM/wind_speed_raster/HadGEM3-GC31-HM/archive.zip
    snakemake -c1 results/input/STORM/RP_raster/HadGEM3-GC31-HM/archive.zip
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/STORM/{EVENTS_OR_RASTERS}/{STORM_MODEL}/archive.zip"
    shell:
        """
        if [[ "{wildcards.EVENTS_OR_RASTERS}" == "events" ]]
        then
            RESOURCE_FILE="events/{wildcards.STORM_MODEL}.txt"
        elif [[ "{wildcards.EVENTS_OR_RASTERS}" == "wind_speed_raster" ]]
        then
            RESOURCE_FILE="wind_speed_raster/{wildcards.STORM_MODEL}.txt"
        else
            RESOURCE_FILE="RP_raster/{wildcards.STORM_MODEL}.txt"
        fi

        wget \
            --input-file=config/hazard_resource_locations/cyclone/STORM/$RESOURCE_FILE \
            --directory-prefix={wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_RASTERS}/{wildcards.STORM_MODEL}/ \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        mv \
            {wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_RASTERS}/{wildcards.STORM_MODEL}/*.zip \
            {wildcards.OUTPUT_DIR}/input/STORM/{wildcards.EVENTS_OR_RASTERS}/{wildcards.STORM_MODEL}/archive.zip
        """

rule parse_storm:
    """
    Process raw CSV track data into common geoparquet format
    Test with:
    snakemake -c1 results/storm_tracks/STORM-constant/0/tracks.geoparquet
    """
    input:
        csv_dir="{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/raw"
    output:
        parquet="{OUTPUT_DIR}/storm_tracks/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet"
    script:
        "./parse_STORM.py"


rule slice_storm:
    """
    Select tracks by location
    To test:
    snakemake -c1 results/power/by_country/PRI/storms/STORM-constant/0/tracks.geoparquet
    """
    input:
        global_tracks="{OUTPUT_DIR}/storm_tracks/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/STORM-{STORM_MODEL}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "./slice_storm_tracks.py"


rule extract_storm_event:
    """
    Unzip a storm file for a basin we are interested in
    Test with:
    snakemake -c1 results/input/STORM/events/constant/NA/STORM_DATA_IBTRACS_NA_1000_YEARS_0.txt
    snakemake -c1 results/input/STORM/events/HadGEM3-GC31-HM/NA/STORM_DATA_HadGEM3-GC31-HM_NA_1000_YEARS_0_IBTRACSDELTA.txt
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


rule extract_all_storm_events:
    """
    Unzip all the storm files for a given model. Rename to appropriate file extension (CSV).
    Test with:
    snakemake -c1 results/input/STORM/events/constant/raw/
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


rule extract_storm_wind_speed_raster:
    """
    Unzip a STORM file of wind speeds (for a given climate model, basin and return period).
    Test with:
    snakemake -c1 results/input/STORM/wind_speed_raster/constant/NA/STORM_FIXED_RETURN_PERIODS_constant_NA_20_YR_RP.tif
    snakemake -c1 results/input/STORM/wind_speed_raster/CMCC-CM2-VHR4/NA/STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_NA_10_YR_RP.tif
    """
    input:
        "{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/archive.zip"
    output:
        "{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    params:
        STORM_MODEL_UPPER = lambda wildcards: wildcards.STORM_MODEL.upper()
    shell:
        """
        unpack_dir="{wildcards.OUTPUT_DIR}/input/STORM/wind_speed_raster/{wildcards.STORM_MODEL}/{wildcards.STORM_BASIN}/"

        if [[ "{wildcards.STORM_MODEL}" == "constant" ]]
        then
            # present climate
            unzip -o {input} STORM_FIXED_RETURN_PERIODS_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif \
                -d {wildcards.OUTPUT_DIR}/input/STORM/wind_speed_raster/constant/{wildcards.STORM_BASIN}/
            mv $unpack_dir/STORM_FIXED_RETURN_PERIODS_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif \
                $unpack_dir/STORM_FIXED_RETURN_PERIODS_constant_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif
        else
            # future climate
            uppercase_filename="STORM_FIXED_RETURN_PERIODS_{params.STORM_MODEL_UPPER}_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif"
            desired_filename="STORM_FIXED_RETURN_PERIODS_{wildcards.STORM_MODEL}_{wildcards.STORM_BASIN}_{wildcards.STORM_RP}_YR_RP.tif"

            # try extracting filename with STORM_MODEL, fall back to STORM_MODEL_UPPER
            unzip -o {input} $desired_filename -d $unpack_dir || unzip -o {input} $uppercase_filename -d $unpack_dir

            # if we unpacked the uppercase version and its different to STORM_MODEL, rename it
            pushd $unpack_dir
                if [ -f $uppercase_filename ] && [ $uppercase_filename != $desired_filename ]; then
                    mv $uppercase_filename $desired_filename
                fi
            popd
        fi
        """


rule extract_storm_RP_raster:
    """
    Unzip a STORM file of return periods in years (for a given climate model, basin and wind speed).
    Test with:
    snakemake -c1 results/input/STORM/RP_raster/constant/NA/STORM_FIXED_WIND_SPEEDS_constant_NA_30_MS.tif
    snakemake -c1 results/input/STORM/RP_raster/CMCC-CM2-VHR4/NA/STORM_FIXED_WIND_SPEEDS_CMCC-CM2-VHR4_NA_30_MS.tif
    snakemake -c1 results/input/STORM/RP_raster/HadGEM3-GC31-HM/NA/STORM_FIXED_WIND_SPEEDS_HadGEM3-GC31-HM_NA_30_MS.tif
    """
    input:
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/archive.zip"
    output:
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_{STORM_BASIN}_{STORM_SPEED}_MS.tif"
    params:
        STORM_MODEL_UPPER = lambda wildcards: wildcards.STORM_MODEL.upper()
    shell:
        """
        unpack_dir="{wildcards.OUTPUT_DIR}/input/STORM/RP_raster/{wildcards.STORM_MODEL}/{wildcards.STORM_BASIN}/"

        if [[ "{wildcards.STORM_MODEL}" == "constant" ]]
        then
            # present climate
            unzip -o {input} STORM_FIXED_WIND_SPEEDS_{wildcards.STORM_BASIN}_{wildcards.STORM_SPEED}_MS.tif \
                -d {wildcards.OUTPUT_DIR}/input/STORM/RP_raster/constant/{wildcards.STORM_BASIN}/
            mv $unpack_dir/STORM_FIXED_WIND_SPEEDS_{wildcards.STORM_BASIN}_{wildcards.STORM_SPEED}_MS.tif \
                $unpack_dir/STORM_FIXED_WIND_SPEEDS_constant_{wildcards.STORM_BASIN}_{wildcards.STORM_SPEED}_MS.tif
        else
            # future climate
            uppercase_filename="STORM_FIXED_WIND_SPEEDS_{params.STORM_MODEL_UPPER}_{wildcards.STORM_BASIN}_{wildcards.STORM_SPEED}_MS.tif"
            desired_filename="STORM_FIXED_WIND_SPEEDS_{wildcards.STORM_MODEL}_{wildcards.STORM_BASIN}_{wildcards.STORM_SPEED}_MS.tif"

            # try extracting filename with STORM_MODEL, fall back to STORM_MODEL_UPPER
            unzip -o {input} $desired_filename -d $unpack_dir || unzip -o {input} $uppercase_filename -d $unpack_dir

            # if we unpacked the uppercase version and its different to STORM_MODEL, rename it
            pushd $unpack_dir
                if [ -f $uppercase_filename ] && [ $uppercase_filename != $desired_filename ]; then
                    mv $uppercase_filename $desired_filename
                fi
            popd
        fi
        """


rule wrap_storm_wind_speed_raster:
    """
    Test with:
    snakemake -c1 results/input/STORM/wind_speed_raster/CMCC-CM2-VHR4/NA/STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_NA_10_YR_RP.wrapped.tif
    """
    input:
        "{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_BASIN}_{STORM_RP}_YR_RP.tif"
    output:
        temp("{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_BASIN}_{STORM_RP}_YR_RP.wrapped.tif")
    shell:
        """
        gdalwarp -te -179.85 -60.15 180.15 59.95 -co COMPRESS=LZW -co PREDICTOR=2 -co TILED=YES {input} {output}
        """


rule wrap_storm_RP_raster:
    """
    Test with:
    snakemake -c1 results/input/STORM/RP_raster/CMCC-CM2-VHR4/NA/STORM_FIXED_WIND_SPEEDS_CMCC-CM2-VHR4_NA_30_MS.wrapped.tif
    """
    input:
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_{STORM_BASIN}_{STORM_SPEED}_MS.tif"
    output:
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/{STORM_BASIN}/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_{STORM_BASIN}_{STORM_SPEED}_MS.wrapped.tif"
    shell:
        """
        gdalwarp -te -179.85 -60.15 180.15 59.95 -co COMPRESS=LZW -co PREDICTOR=2 -co TILED=YES -dstnodata 0 {input} {output}
        """


rule mosaic_storm_wind_speed_raster:
    """
    Merge basin return period maps to global extent
    Test with:
    snakemake -c1 results/input/STORM/wind_speed_raster/constant/STORM_FIXED_RETURN_PERIODS_constant_10_YR_RP.tif
    snakemake -c1 results/input/STORM/wind_speed_raster/CMCC-CM2-VHR4/STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_10_YR_RP.tif
    """
    input:
        basin_tif_ep="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/EP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_EP_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_na="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/NA/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_NA_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_ni="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/NI/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_NI_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_si="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/SI/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_SI_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_sp="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/SP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_SP_{STORM_RP}_YR_RP.wrapped.tif",
        basin_tif_wp="{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/WP/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_WP_{STORM_RP}_YR_RP.wrapped.tif",
    output:
        "{OUTPUT_DIR}/input/STORM/wind_speed_raster/{STORM_MODEL}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_RP}_YR_RP.tif"
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
            --NoDataValue=0 \
            --creation-option="COMPRESS=LZW" \
            --creation-option="PREDICTOR=2" \
            --creation-option="TILED=YES"
        """


rule mosaic_storm_RP_raster:
    """
    Merge basin wind speed maps to global extent
    Test with:
    snakemake -c1 results/input/STORM/RP_raster/constant/STORM_FIXED_WIND_SPEEDS_constant_30_MS.tif
    snakemake -c1 results/input/STORM/RP_raster/CMCC-CM2-VHR4/STORM_FIXED_WIND_SPEEDS_CMCC-CM2-VHR4_30_MS.tif
    """
    input:
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/EP/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_EP_{STORM_SPEED}_MS.wrapped.tif",
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/NA/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_NA_{STORM_SPEED}_MS.wrapped.tif",
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/NI/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_NI_{STORM_SPEED}_MS.wrapped.tif",
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/SI/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_SI_{STORM_SPEED}_MS.wrapped.tif",
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/SP/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_SP_{STORM_SPEED}_MS.wrapped.tif",
        "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/WP/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_WP_{STORM_SPEED}_MS.wrapped.tif",
    output:
        mosaic = "{OUTPUT_DIR}/input/STORM/RP_raster/{STORM_MODEL}/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_{STORM_SPEED}_MS.tif"
    run:
        import numpy as np
        import rioxarray as rio
        import xarray as xr

        max_return_period_years = 2_000  # discard values above this threshold

        rasters: list[xr.DataArray] = []
        for basin_path in input:
            basin: str = basin_path.split("/")[-2]
            raster: xr.DataArray = rio.open_rasterio(basin_path).sel(band=1)
            raster.expand_dims("basin")
            raster["basin"] = basin
            # discard return period values over some threshold (those with only a handful of events)
            raster = raster.where(raster < max_return_period_years)
            # we're interested in the minimum non-zero value across the (overlapping) basins
            # set zero values to NaN, allows nanmin subsequently
            raster = raster.where(raster != 0)
            rasters.append(raster)
        concat: xr.DataArray = xr.concat(rasters, "basin")

        mosaic: xr.DataArray = raster.copy(data=np.nanmin(concat.values, axis=0))
        mosaic = mosaic.drop_vars("basin")
        mosaic.rio.to_raster(output.mosaic)


rule mosaic_storm_wind_speed_raster_all:
    """
    Target rule to create mosaicked wind speed rasters for all models and return periods
    Test with:
    snakemake -c1 results/input/STORM/wind_speed_raster/mosaic.done
    """
    input:
        tiffs=expand(
            "{{OUTPUT_DIR}}/input/STORM/wind_speed_raster/{STORM_MODEL}/STORM_FIXED_RETURN_PERIODS_{STORM_MODEL}_{STORM_RP}_YR_RP.tif",
            STORM_MODEL=[
                "constant",
                "CMCC-CM2-VHR4",
                "CNRM-CM6-1-HR",
                "EC-Earth3P-HR",
                "HadGEM3-GC31-HM",
            ],
            STORM_RP=(
                list(range(10, 100, 10))
                + list(range(100, 1000, 100))
                + list(range(1000, 10001, 1000))
            ),
        )
    output:
        touch("{OUTPUT_DIR}/input/STORM/wind_speed_raster/mosaic.done")


rule mosaic_storm_RP_raster_all:
    """
    Target rule to create mosaicked return period rasters for all models and return periods
    N.B. Rule mosaic_storm_RP_raster drops return period values over a threshold (considered noise).
    Test with:
    snakemake -c1 results/input/STORM/RP_raster/mosaic.done
    """
    input:
        tiffs=expand(
            "{{OUTPUT_DIR}}/input/STORM/RP_raster/{STORM_MODEL}/STORM_FIXED_WIND_SPEEDS_{STORM_MODEL}_{STORM_SPEED}_MS.tif",
            STORM_MODEL=[
                "constant",
                "CMCC-CM2-VHR4",
                "CNRM-CM6-1-HR",
                "EC-Earth3P-HR",
                "HadGEM3-GC31-HM",
            ],
            STORM_SPEED=(20, 25, 29, 30, 35, 37, 40, 43, 45, 50, 51, 55, 60, 61, 65, 70, 75),
        )
    output:
        touch("{OUTPUT_DIR}/input/STORM/RP_raster/mosaic.done")
