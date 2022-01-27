"""Extracts winds data for specified region sample and year.

"""
import os


def required_years_remaining(YRS):
    years_completed_files = glob(os.path.join(DATA_DIR, "intersection", "storm_data", "all_winds", "*csv"))
    years_completed = [int(file[file.find("_y")+2:-4]) for file in years_completed_files]  # single out years from directory files
    years_remaining = [int(YR) for YR in YRS if int(YR) not in years_completed]  # list all years that are not already completed (years_completed)
    print('Years to cover: ',years_remaining)
    return years_remaining

YEARS_req = required_years_remaining(YEARS)

TC_years = expand(
    os.path.join(
        "data",
        "intersection",
        "storm_data",
        "all_winds",
        "log",
        "__winds_completed_r{region}_s{sample}_y{year}.txt",
    ),
    region=REGIONS,
    sample=SAMPLES,
    year=YEARS_req,
)


rule intersect_winds_indiv:
    input:
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
        os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
        os.path.join(
            "data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"
        ),
        os.path.join(
            "data",
            "stormtracks",
            "events",
            "STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
        ),
    output:
        [
            os.path.join(
                "data",
                "intersection",
                "storm_data",
                "all_winds",
                "log",
                "__winds_completed_r{region}_s{sample}" + f"_y{year}.txt",
            )
            for year in YEARS_req
        ],  # use as pseudo file for csv (otherwise snakemake deletes output files everytime)
    shell:
        (
            "python3 "
            + os.path.join("workflow", "scripts", "intersect", "intersect_3_winds.py")
            + " {wildcards.region} {wildcards.sample} "
            + '"""'
            + str(YEARS)
            + '"""'
        )


rule intersect_wind:
    input:
        TC_years,
