"""
The TrackSet container: track points plus the metadata needed to use them.

A bag of track points is not enough to do risk analysis with: you also need
to know the time span the set represents (to turn event counts into
frequencies), whether its timestamps are real or fabricated, and what wind
speed averaging period it reports. Upstream formats carry this information
out-of-band (in filenames, READMEs and folklore); TrackSet carries it with
the data, including through GeoParquet serialisation.
"""

import dataclasses
import io
import json
from dataclasses import dataclass, field
from os import PathLike
from typing import Any, Union

import geopandas as gpd
import pyarrow.parquet as pq
from shapely.geometry.base import BaseGeometry

from . import ops, schema

_PARQUET_METADATA_KEY = b"trackset"


@dataclass(frozen=True)
class TrackSet:
    """
    A set of tropical cyclone tracks in the canonical schema, with metadata.
    """

    #: Track points conforming to :mod:`trackset.schema`.
    data: gpd.GeoDataFrame
    #: Track set family and variant, e.g. "IBTrACS", "STORM-constant",
    #: "CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050".
    source: str
    #: Total years of observation or simulation this set represents. The
    #: denominator for annual frequencies -- NOT necessarily the span of the
    #: ``year`` column (independent samples of the same epoch accumulate).
    years: float
    #: True if the DatetimeIndex is fabricated (synthetic sets report only
    #: timestep numbers; see :func:`trackset.ops.synthetic_time_index`).
    synthetic_time: bool
    #: Averaging period of ``max_wind_speed_ms``. Readers normalise to 1-minute
    #: sustained winds on ingest.
    wind_averaging_period: str = "1min"
    #: Free-form provenance: scenario, GCM, epoch, upstream file names, etc.
    attributes: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        schema.validate(self.data)
        if self.years <= 0:
            raise ValueError(f"years must be positive, got {self.years}")

    @property
    def n_tracks(self) -> int:
        return self.data["track_id"].nunique()

    @property
    def annual_frequency(self) -> float:
        """Mean number of tracks per year represented by this set."""
        return self.n_tracks / self.years

    def subset(self, geometry: BaseGeometry, buffer_deg: float = 0.0) -> "TrackSet":
        """
        Subset to tracks passing within ``buffer_deg`` of ``geometry``. See
        :func:`trackset.ops.subset_by_geometry`. Metadata (including
        ``years``: the unselected tracks still happened) is preserved.
        """

        return dataclasses.replace(
            self, data=ops.subset_by_geometry(self.data, geometry, buffer_deg)
        )

    def to_parquet(self, path: Union[str, PathLike]) -> None:
        """
        Write to GeoParquet with TrackSet metadata embedded in the file-level
        key-value metadata (under the ``trackset`` key).
        """

        buffer = io.BytesIO()
        self.data.to_parquet(buffer)
        buffer.seek(0)
        table = pq.read_table(buffer)

        metadata = dict(table.schema.metadata or {})
        metadata[_PARQUET_METADATA_KEY] = json.dumps(
            {
                "source": self.source,
                "years": self.years,
                "synthetic_time": self.synthetic_time,
                "wind_averaging_period": self.wind_averaging_period,
                "attributes": self.attributes,
            }
        ).encode()
        pq.write_table(table.replace_schema_metadata(metadata), path)

    @classmethod
    def read_parquet(cls, path: Union[str, PathLike]) -> "TrackSet":
        """
        Read a GeoParquet file written by :meth:`to_parquet`.
        """

        file_metadata = pq.read_schema(path).metadata or {}
        if _PARQUET_METADATA_KEY not in file_metadata:
            raise ValueError(
                f"{path} has no '{_PARQUET_METADATA_KEY.decode()}' metadata; "
                "was it written by TrackSet.to_parquet?"
            )
        metadata = json.loads(file_metadata[_PARQUET_METADATA_KEY])

        return cls(data=gpd.read_parquet(path), **metadata)
