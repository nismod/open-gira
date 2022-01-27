"""Merges the gdp loss statistics.

"""
import os

stat_csv = os.path.join("data", "intersection", "combined_storm_statistics.csv")


rule merge_all_stats:
    """Use this rule for the combined stats file"""
    input:
        TC_years,
        region_grid,
        expand(
            os.path.join(
                "data",
                "intersection",
                "storm_data",
                "damages",
                "storm_r{region}_s{sample}_y{year}.txt",
            ),
            region=REGIONS,
            sample=SAMPLES,
            year=YEARS,
        ),
    output:
        stat_csv,
    shell:
        "python3 " + os.path.join("workflow", "scripts", "intersect", "stat_merger.py")
