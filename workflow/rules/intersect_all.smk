"""Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""


rule intersect_all:
    input:
        TC_all,
        stat_csv,
        storm_details_all