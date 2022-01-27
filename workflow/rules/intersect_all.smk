"""Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""


rule intersect_all:
    input:
        TC_years,  # required for stat_csv to be up to date
        region_grid,  # required for stat_csv to be up to date
        stat_csv,
