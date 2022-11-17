"""
Field names used throughout data pipelines.
"""

# exposure table hazard intensity fields expected to be prefixed as such
HAZARD_PREFIX = "hazard-"

# exposure table field containing the cost to rebuild per unit length
REHAB_COST = "rehab_cost_USD_per_km"

# length of edge calculated after intersection
SPLIT_LENGTH = "length_km"

# indicies of relevant cell in associated raster file(s)
RASTER_I = "raster_i"
RASTER_J = "raster_j"
