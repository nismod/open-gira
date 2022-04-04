"""Merges the gdp loss statistics.

"""
import os

stat_csv = os.path.join("data", "intersection", "combined_storm_statistics.csv")


storm_details_all = expand(
    [
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "individual_storms",
            "{region}",
            "storm_{nh}",
            "storm_r{region}_s{sample}_n{nh}.txt",
        ),
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "individual_storms",
            "{region}",
            "storm_{nh}",
            "storm_track_r{region}_s{sample}_n{nh}.gpkg",
        ),
    ],
    region=REGIONS,
    sample=SAMPLES,
    nh=find_nh_mult(YEARS, REGIONS, SAMPLES)[2],
)
#print(storm_details_all)

rule merge_all_stats:
    """Use this rule for the combined stats file"""
    input:
        storm_details_all,
    output:
        stat_csv,
    script:
            os.path.join("..", "..", "scripts", "intersect", "stat_merger.py"
        )