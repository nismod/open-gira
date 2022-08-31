"""
Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""

completed_files = expand(
    os.path.join(
        config["output_dir"],
        "power_intersection",
        "storm_data",
        "all_winds",
        "{region}",
        "{sample}",
        "{region}_{sample}_completed.txt",
    ),
    sample=SAMPLES,
    region=REGIONS,
)


rule intersect_all:
    input:
        STORM_STATS_BY_THRESHOLD,
        STORM_STATS_BY_REGION_SAMPLE_THRESHOLD,  # to ensure all stat files up to date and present (not just the final csv (above))
        completed_files,  # to ensure all storm files up to date and present
