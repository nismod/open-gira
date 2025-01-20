import geopandas
import numpy as np
import rasterio.features


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


#
# Raster manipulation helpers (working with 2d numpy arrays)
#
def clip_array(arr, block_size):
    """Clip a 2d array to an integer multiple of block_size in each
    dimension"""
    clip_rows = arr.shape[0] - (arr.shape[0] % block_size)
    clip_cols = arr.shape[1] - (arr.shape[1] % block_size)

    clipped = arr[0:clip_rows, 0:clip_cols]
    return clipped


def resample_sum(arr, block_size):
    """Resample a 2d array, summing each block of (block_size x
    block_size) to give each cell in the output array"""
    nblocks_0 = arr.shape[0] // block_size
    nblocks_1 = arr.shape[1] // block_size

    blocks = arr.reshape(nblocks_0, block_size, nblocks_1, block_size)

    return np.sum(blocks, axis=(1, 3))


def repeat_2d(arr, block_size):
    """Repeat each element in a 2d array, so each value fills a (block_size x
    block_size) area"""
    return np.repeat(np.repeat(arr, block_size, axis=0), block_size, axis=1)


def floor_int(a):
    """Floor and convert to integer"""
    return np.floor(a).astype(int)


def zero_divide(a, b):
    """Divide (a / b) but return zero where (b == 0)"""
    return np.divide(a, b, out=np.zeros_like(a, dtype="float64"), where=(b != 0))


def rasterize(gdf: geopandas.GeoDataFrame, column: str, template_ds):
    """Burn values from a GeoDataFrame column into a raster of shape and transform
    specified by template_ds"""
    return rasterio.features.rasterize(
        (
            (f["geometry"], f["properties"][column])
            for f in gdf.__geo_interface__["features"]
        ),
        out_shape=template_ds.shape,
        transform=template_ds.transform,
    )
