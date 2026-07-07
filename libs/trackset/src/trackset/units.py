"""
Unit conversion constants and wind speed averaging-period adjustments.
"""

MS_PER_KNOT = 0.5144
KM_PER_NAUTICAL_MILE = 1.852

# Divide by this factor to 'convert' 10-minutely sustained winds to 1-minutely
# sustained wind speeds, noting the vagaries of this process as explained here:
# https://library.wmo.int/doc_num.php?explnum_id=290
TEN_MINUTE_TO_ONE_MINUTE_WIND_FACTOR = 0.88

# Wind scale and shift values taken from CLIMADA: climada.hazard.tc_tracks
#
# Used to 'convert' between the reporting agencies' averaging periods and a
# 1-minutely period: subtract the shift, divide by the scale.
#
# From Table 1 in: Knapp, K.R., Kruk, M.C. (2010): Quantifying Interagency
# Differences in Tropical Cyclone Best-Track Wind Speed Estimates.
# Monthly Weather Review 138(4): 1459-1473.
# https://library.wmo.int/index.php?lvl=notice_display&id=135
IBTRACS_AGENCY_1MIN_WIND_FACTOR: dict[str, tuple[float, float]] = {
    "USA": (1.0, 0.0),
    "TOKYO": (0.60, 23.3),
    "NEWDELHI": (1.0, 0.0),
    "REUNION": (0.88, 0.0),
    "BOM": (0.88, 0.0),
    "NADI": (0.88, 0.0),
    "WELLINGTON": (0.88, 0.0),
    "CMA": (0.871, 0.0),
    "HKO": (0.9, 0.0),
    "DS824": (1.0, 0.0),
    "TD9636": (1.0, 0.0),
    "TD9635": (1.0, 0.0),
    "NEUMANN": (0.88, 0.0),
    "MLC": (1.0, 0.0),
}
