"""
Exporters to other tropical cyclone risk tools.
"""

from .core import TrackSet


def to_climada(trackset: TrackSet):
    """
    Convert to a ``climada.hazard.TCTracks`` instance.

    Not yet implemented. The intended mapping, per track:

    - one ``xarray.Dataset`` per ``track_id``, with ``time`` coordinate from
      the index and ``lat``/``lon`` coordinates from the geometry
    - ``max_sustained_wind`` from ``max_wind_speed_ms`` (CLIMADA convention:
      1-minute sustained, unit recorded in attrs)
    - ``central_pressure`` from ``min_pressure_hpa``
    - ``radius_max_wind`` from ``radius_to_max_winds_km`` (converted to
      nautical miles, CLIMADA's convention)
    - ``basin``, ``category``, ``name`` attributes where present
    - ``TrackSet.years`` informing per-track ``frequency`` (CLIMADA's annual
      occurrence rate: 1 / years for each simulated track)

    climada is an optional dependency: install ``climada`` to use this once
    implemented.
    """

    raise NotImplementedError("planned: see docstring for the intended field mapping")
