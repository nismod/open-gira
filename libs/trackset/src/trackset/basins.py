"""
Tropical cyclone basin definitions and point tagging.

Basin extents follow the STORM dataset definitions (Bloemendaal et al.):
https://www.nature.com/articles/s41597-020-0381-2/tables/2

The polygons are defined in 0-360 longitude space so that no basin crosses
the antimeridian (the South Pacific does in signed coordinates);
:func:`tag_basin` handles the shift for signed-longitude points.
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

_DELTA = 1e-3

#: Basin boundary vertices in (longitude, latitude), signed degrees.
_BASIN_BOUNDARIES: tuple[tuple[str, tuple[tuple[float, float], ...]], ...] = (
    (
        "EP",
        (
            (-180, 60),
            (-180, 0),
            (-75, 0),
            (-75, 10),
            (-85, 10),
            (-85, 15),
            (-90, 15),
            (-90, 17.5),
            (-100, 17.5),
            (-100, 60),
        ),
    ),
    (
        "NA",
        (
            (-100, 60),
            (-100, 17.5),
            (-90, 17.5),
            (-90, 15),
            (-85, 15),
            (-85, 10),
            (-75, 10),
            (-75, 0),
            (360 - _DELTA, 0),
            (360 - _DELTA, 60),
        ),
    ),
    ("NI", ((30, 60), (30, 0), (100, 0), (100, 60))),
    ("SI", ((10, 0), (10, -60), (135, -60), (135, 0))),
    ("SP", ((135, 0), (135, -60), (-120, -60), (-120, 0))),
    ("WP", ((100, 60), (100, 0), (180, 0), (180, 60))),
)


def _signed_longitude_to_strictly_positive(coords):
    return [(lon + 360 if lon < 0 else lon, lat) for lon, lat in coords]


def basin_polygons() -> gpd.GeoDataFrame:
    """
    The tropical cyclone basins as polygons in 0-360 longitude space, with a
    ``basin_id`` column ("EP", "NA", "NI", "SI", "SP", "WP").
    """

    basin_ids, boundaries = zip(*_BASIN_BOUNDARIES)
    return gpd.GeoDataFrame(
        data={
            "basin_id": list(basin_ids),
            "geometry": [
                Polygon(_signed_longitude_to_strictly_positive(boundary))
                for boundary in boundaries
            ],
        },
        crs=4326,
    )


def tag_basin(df, lon, lat, basins: gpd.GeoDataFrame | None = None):
    """
    Return a copy of ``df`` with a ``basin_id`` column from a spatial join of
    point coordinates against the basin polygons.

    Longitudes may be signed or 0-360. Points outside every basin (beyond 60
    degrees of latitude) are dropped, matching the upstream parsers.

    Args:
        df: Any DataFrame of track points.
        lon: Point longitudes, degrees.
        lat: Point latitudes, degrees.
        basins: Basin polygons; defaults to :func:`basin_polygons`.
    """

    if basins is None:
        basins = basin_polygons()

    lon = np.asarray(lon, dtype=float)
    lon = np.where(lon < 0, lon + 360, lon)
    points = gpd.GeoDataFrame(
        df.copy(),
        geometry=gpd.points_from_xy(lon, lat),
        crs=4326,
    )
    tagged = points.sjoin(basins.to_crs(epsg=4326)).drop(columns="index_right")
    return tagged.drop(columns="geometry")
