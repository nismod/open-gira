"""Download STORM IBTrACS present climate synthetic tropical cyclone tracks and tropical cyclone 
wind speed return periods

Reference
---------
https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085
https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164
"""

out_fixed = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "stormtracks",
        "fixed",
        "STORM_FIXED_{param}_{region}.nc",
    ),
    region=config["regions"],
    param=["RETURN_PERIODS", "TC_WIND_SPEEDS"],
)

storm_model_type = config["storm_model_type"]

if storm_model_type == "constant":
    wind_file_start = "STORM_DATA_IBTRACS_"
    wind_file_end = ""
    unzip_file = "STORM_DATA3.zip"
elif storm_model_type == "CMCC-CM2-VHR4":
    wind_file_start = "STORM_DATA_CMCC-CM2-VHR4_"
    wind_file_end = "_IBTRACSDELTA"
    unzip_file = "CMCC"
elif storm_model_type == "CNRM-CM6-1-HR":
    wind_file_start = "STORM_DATA_CNRM-CM6-1-HR_"
    wind_file_end = "_IBTRACSDELTA"
    unzip_file = "CNRM"
elif storm_model_type == "EC-Earth3P-HR":
    wind_file_start = "STORM_DATA_EC-Earth3P-HR_"
    wind_file_end = "_IBTRACSDELTA"

    unzip_file = "ECEARTH"
elif storm_model_type == "HadGEM3-GC31-HM":
    wind_file_start = "STORM_DATA_HadGEM3-GC31-HM_"
    wind_file_end = "_IBTRACSDELTA"
    unzip_file = "HADGEM"
else:
    raise RuntimeError(
        f"The selected storm type model ({storm_model_type}) is not a valid option"
    )

out_events = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "stormtracks",
        "events",
        wind_file_start + "{region}_1000_YEARS_{num}" + wind_file_end + ".txt",
    ),
    region=CYCLONE_REGIONS,
    num=list(range(0, 10)),
)


rule download_stormtracks_fixed:
    output:
        out_fixed,
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
        out_events,
    shell:
        f"""
        wget \
            --input-file=workflow/scripts/storm_tracks_{storm_model_type}.txt \
            --directory-prefix={config['output_dir']}/input/stormtracks/events \
            --timestamping \
            --no-check-certificate \
            --content-disposition
        unzip -o {config['output_dir']}/input/stormtracks/events/{unzip_file} -d {config['output_dir']}/input/stormtracks/events
        """
