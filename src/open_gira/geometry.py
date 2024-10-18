import geopandas


def buffer_via_reprojection(
    geoms: geopandas.GeoSeries, buffer_radius_m
) -> geopandas.GeoSeries:
    """Buffer geographical geometries

    First project into a UTM CRS, estimated based on the bounds of the dataset,
    then buffer, then project back into original CRS.

    Parameters
    ----------
    geoms: GeoDataFrame
    buffer_radius_m: float
    """
    projected_crs = geoms.estimate_utm_crs()
    return (
        geoms.geometry.to_crs(projected_crs).buffer(buffer_radius_m).to_crs(geoms.crs)
    )
