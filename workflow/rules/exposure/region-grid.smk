"""Intersect power network with return period maps

- split edges on global grid
- assign edges and nodes to grid indices
- join windspeed per return period map

Next:
- calculate winds per event, only needed for unique grid indices with infrastructure
"""

import os


rule intersect_unit_generator:
    conda: "../../../environment.yml"
    input:
        network="{OUTPUT_DIR}/processed/power/{BOX}/edges_{BOX}.parquet",
        tifs="{OUTPUT_DIR}/input/storm-ibtracs/fixed/*"
    output:
        geoparquet="{OUTPUT_DIR}/processed/power/{BOX}/edges_split_{BOX}.geoparquet",
    script:
        "../../scripts/intersection.py"
