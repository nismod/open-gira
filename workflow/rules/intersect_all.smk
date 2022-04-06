"""Perform intersection

Intersect the infrastructure data with the hazard (storms)
"""


rule intersect_all:
    input:
        stat_csv,
