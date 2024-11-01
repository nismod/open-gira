import geopandas
import numpy
import pandas
import rasterstats


def max_vector_rasters_intersection(
    vector: geopandas.GeoDataFrame, rasters: pandas.DataFrame
) -> geopandas.GeoDataFrame:
    """Intersect vector geometries with raster files, adding columns with
    maximum values from the intersection

    Parameters
    ----------
    vector: GeoDataFrame
        vector geometries. Output columns will be added to this GeoDataFrame.
    rasters: DataFrame
        metadata table with "key" and "path" columns. "key" is used for output
        column names. "path" is used to specify raster file paths.
    """
    for raster in rasters.itertuples():
        vector[raster.key] = max_vector_raster_intersection(
            vector.geometry, raster.path
        )
    return vector


def max_vector_raster_intersection(
    vector: geopandas.GeoSeries, raster: str
) -> numpy.array:
    """Intersect vector geometries with raster, return array of max raster values

    Parameters
    ----------
    vector: GeoSeries
        vector (point/line/polygon) geometries
    raster: str
        path to raster file
    """
    maxs = numpy.zeros(len(vector.geometry))
    for i, stats in enumerate(
        rasterstats.gen_zonal_stats(vector, raster, stats="max", all_touched=True)
    ):
        maxs[i] = stats["max"]
    return maxs
