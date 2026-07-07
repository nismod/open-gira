"""
trackset: tropical cyclone track sets, normalised.

Read observed (IBTrACS) and synthetic (STORM, IRIS, CHAZ, Emanuel) tropical
cyclone track sets into a single, validated GeoDataFrame schema, with the
metadata needed for risk analysis (time span represented, wind averaging
period, provenance) carried alongside the data and through GeoParquet.
"""

from . import ops, units
from .core import TrackSet
from .readers import (
    read_chaz,
    read_emanuel,
    read_ibtracs,
    read_iris,
    read_storm,
)
from .schema import OPTIONAL, REQUIRED, SchemaError, problems, validate

__version__ = "0.1.0a1"

__all__ = [
    "TrackSet",
    "SchemaError",
    "REQUIRED",
    "OPTIONAL",
    "validate",
    "problems",
    "ops",
    "units",
    "read_chaz",
    "read_emanuel",
    "read_ibtracs",
    "read_iris",
    "read_storm",
]
