"""
Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""

rule intersect_all:
    input:
        STORM_STATS_BY_THRESHOLD,
        STORM_STATS_BY_REGION_SAMPLE_THRESHOLD,  # to ensure all stat files up to date and present (not just the final csv (above))
        COMPLETION_FLAG_FILES,  # to ensure all storm files up to date and present
