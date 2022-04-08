"""Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""


rule intersect_all:
    input:
        stat_csv,
        all_indiv_stat_csv  # to ensure all files up to date and present (not just the final csv (above))
