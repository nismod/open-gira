import numpy as np
import rasterio

from affine import Affine
from snail.damages import PiecewiseLinearDamageCurve


def clip_array(arr, block_size):
    clip_rows = arr.shape[0] - (arr.shape[0] % block_size)
    clip_cols = arr.shape[1] - (arr.shape[1] % block_size)

    clipped = arr[0:clip_rows, 0:clip_cols]
    return clipped


def resample_sum(arr, block_size):
    nblocks_0 = arr.shape[0] // block_size
    nblocks_1 = arr.shape[1] // block_size

    blocks = arr.reshape(nblocks_0, block_size, nblocks_1, block_size)

    return np.sum(blocks, axis=(1, 3))


def repeat_2d(arr, block_size):
    """Repeat each element from a 2d array, so each value fills a (block_size x block_size) area"""
    return np.repeat(np.repeat(arr, block_size, axis=0), block_size, axis=1)


def read_ds(ds, band=1, replace_nodata=False, nodata_fill=0):
    data = ds.read(band)
    if replace_nodata:
        data = np.where(data == ds.nodata, nodata_fill, data)
    return data


def to_int(a):
    return np.floor(a).astype(int)


def main(value_150ss_tif, volume_3ss_tif, flood_1ss_tif, prefix):
    with rasterio.open(value_150ss_tif) as value_150ss_ds:
        value_150ss_all = read_ds(value_150ss_ds, replace_nodata=True)

    with rasterio.open(volume_3ss_tif) as volume_3ss_ds:
        volume_3ss_all = read_ds(volume_3ss_ds, replace_nodata=True)

    # lon, lat of volume_3ss top left
    volume_3ss_all_ul_xy = volume_3ss_ds.transform * (0, 0)
    # col, row in value_150ss_all, inset one extra
    value_150ss_ul_cr = to_int(~value_150ss_ds.transform * (volume_3ss_all_ul_xy)) + 1
    # lon, lat of that value_150ss_all pixel - this is our new top left
    ul_xy_150ss = value_150ss_ds.transform * value_150ss_ul_cr
    # col, row in volume_3ss_all
    volume_3ss_ul_cr = to_int(~volume_3ss_ds.transform * ul_xy_150ss)
    # lon, lat of that volume_3ss_all pixel - new top left for 3ss purposes (tiny bit offset)
    ul_xy_3ss = volume_3ss_ds.transform * volume_3ss_ul_cr
    ul_xy_150ss, ul_xy_3ss

    # Clip out volume array
    col_idx, row_idx = volume_3ss_ul_cr
    volume_3ss = volume_3ss_all[row_idx:, col_idx:]
    volume_3ss = clip_array(volume_3ss, 50)
    # Resample volume to coarse-scale, "sum"
    volume_150ss = resample_sum(volume_3ss, 50)
    volume_150ss.shape

    # Adapt transform to new top-left and resolution
    a, b, c, d, e, f = volume_3ss_ds.transform[:6]
    t_150ss = Affine(a * 50, b, ul_xy_150ss[0], d, e * 50, ul_xy_150ss[1])
    t_3ss = Affine(a, b, ul_xy_3ss[0], d, e, ul_xy_3ss[1])
    t_150ss, t_3ss

    col_idx, row_idx = value_150ss_ul_cr
    ncols, nrows = volume_150ss.shape
    value_150ss = value_150ss_all[col_idx : col_idx + ncols, row_idx : row_idx + nrows]

    with rasterio.open(
        f"input/giri/THA/{prefix}_vol_150ss.tif",
        "w",
        driver="GTiff",
        height=volume_150ss.shape[0],
        width=volume_150ss.shape[1],
        count=1,
        dtype="float64",
        crs="+proj=latlong",
        transform=t_150ss,
    ) as ds:
        ds.write(volume_150ss, indexes=1)

    with rasterio.open(
        f"input/giri/THA/{prefix}_vol_3ss.tif",
        "w",
        driver="GTiff",
        height=volume_3ss.shape[0],
        width=volume_3ss.shape[1],
        count=1,
        dtype=volume_3ss.dtype,
        crs="+proj=latlong",
        transform=t_3ss,
    ) as ds:
        ds.write(volume_3ss, indexes=1)

    if value_150ss.shape != volume_150ss.shape:
        print("CHKS", value_150ss.shape, volume_150ss.shape)
        assert False

    # Calculate value per unit volume
    # value_per_volume_150ss = value_150ss / volume_150ss
    value_per_volume_150ss = np.divide(
        value_150ss,
        volume_150ss,
        out=np.zeros_like(value_150ss),
        where=volume_150ss != 0,
    )
    # Resample to fine-scale value per volume, "nearest"
    value_per_volume_3ss = repeat_2d(value_per_volume_150ss, 50)
    # Calculate fine-scale value
    value_3ss = value_per_volume_3ss * volume_3ss

    with rasterio.open(
        f"input/giri/THA/{prefix}_val_vol_150ss.tif",
        "w",
        driver="GTiff",
        height=value_per_volume_150ss.shape[0],
        width=value_per_volume_150ss.shape[1],
        count=1,
        dtype=value_per_volume_150ss.dtype,
        crs="+proj=latlong",
        transform=t_150ss,
    ) as ds:
        # Write to window
        ds.write(value_per_volume_150ss, indexes=1)

    with rasterio.open(
        f"input/giri/THA/{prefix}_val_vol_3ss.tif",
        "w",
        driver="GTiff",
        height=value_per_volume_3ss.shape[0],
        width=value_per_volume_3ss.shape[1],
        count=1,
        dtype=value_per_volume_3ss.dtype,
        crs="+proj=latlong",
        transform=t_3ss,
    ) as ds:
        # Write to window
        ds.write(value_per_volume_3ss, indexes=1)

    with rasterio.open(
        f"input/giri/THA/{prefix}_val_3ss.tif",
        "w",
        driver="GTiff",
        height=value_3ss.shape[0],
        width=value_3ss.shape[1],
        count=1,
        dtype=value_3ss.dtype,
        crs="+proj=latlong",
        transform=t_3ss,
    ) as ds:
        # Write to window
        ds.write(value_3ss, indexes=1)

    #
    # Flood intersection
    #
    with rasterio.open(flood_1ss_tif, "r") as flood_1ss_ds:
        flood_1ss = read_ds(flood_1ss_ds, replace_nodata=True)

    # lon, lat of footprint top left
    flood_1ss_ul_xy = flood_1ss_ds.transform * (0, 0)
    # col, row in value_3ss
    t_3ss_ul_cr = to_int(~t_3ss * (flood_1ss_ul_xy))
    # lon, lat of that pixel - this is our new top left
    footprint_ul_xy_3ss = t_3ss * t_3ss_ul_cr
    # col, row in flood_1ss
    flood_1ss_ul_cr = to_int(~flood_1ss_ds.transform * footprint_ul_xy_3ss)
    # lon, lat of that flood_1ss pixel - new top left for 1ss purposes (tiny bit offset)
    ul_xy_1ss = flood_1ss_ds.transform * flood_1ss_ul_cr
    flood_1ss_ul_xy, footprint_ul_xy_3ss, ul_xy_1ss

    # TODO should new top left be greater, not less, in both x and y values?

    # clip to match coarser array extent
    flood_1ss_clipped = clip_array(flood_1ss, 3)
    flood_1ss_height, flood_1ss_width = flood_1ss_clipped.shape

    # lon, lat of footprint lower right
    flood_1ss_lr_xy = flood_1ss_ds.transform * (flood_1ss_width, flood_1ss_height)
    # col, row in value_3ss
    t_3ss_lr_cr = to_int(~t_3ss * (flood_1ss_lr_xy))

    ulc, ulr = t_3ss_ul_cr
    lrc, lrr = t_3ss_lr_cr
    footprint_value_3ss = value_3ss[ulr:lrr, ulc:lrc]

    footprint_value_1ss = repeat_2d(footprint_value_3ss, 3) / 9

    if prefix == "res":
        curve_file = "../bundled_data/damage_curves/flood/residential_asia.csv"
    else:
        curve_file = "../bundled_data/damage_curves/flood/commercial_asia.csv"

    building_flood_depth_damage_curve = PiecewiseLinearDamageCurve.from_csv(
        curve_file,
        intensity_col="inundation_depth_(m)",
        damage_col="damage_fraction",
    )

    if footprint_value_1ss.shape != flood_1ss_clipped.shape:
        print("CHKS", footprint_value_1ss.shape, flood_1ss_clipped.shape)
        assert False

    damage_fraction_1ss = building_flood_depth_damage_curve.damage_fraction(
        flood_1ss_clipped
    )
    damage_value_1ss = footprint_value_1ss * damage_fraction_1ss

    # Adapt transform to new top-left and resolution
    a, b, c, d, e, f = flood_1ss_ds.transform[:6]
    t_1ss = Affine(a, b, ul_xy_1ss[0], d, e, ul_xy_1ss[1])
    t_1ss

    with rasterio.open(
        f"input/giri/THA/{prefix}_dmg_frac_1ss.tif",
        "w",
        driver="GTiff",
        height=damage_fraction_1ss.shape[0],
        width=damage_fraction_1ss.shape[1],
        count=1,
        dtype=damage_fraction_1ss.dtype,
        crs="+proj=latlong",
        transform=t_1ss,
    ) as ds:
        ds.write(damage_fraction_1ss, indexes=1)

    with rasterio.open(
        f"input/giri/THA/{prefix}_dmg_val_1ss.tif",
        "w",
        driver="GTiff",
        height=damage_value_1ss.shape[0],
        width=damage_value_1ss.shape[1],
        count=1,
        dtype=damage_value_1ss.dtype,
        crs="+proj=latlong",
        transform=t_1ss,
    ) as ds:
        ds.write(damage_value_1ss, indexes=1)

    """
    ADM1 damage values:

        exactextract \
            -p ../../admin-boundaries/tha_adm1.shp \
            -r res_dmg_val_1ss.tif \
            -f GID_1 \
            -s sum \
            -o res_dmg_val_1ss.csv

    ADM1 total built volume:

        exactextract \
            -p ../../admin-boundaries/tha_adm1.shp \
            -r ../../ghsl/THA/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0__THA.tif \
            -f GID_1 \
            -s sum \
            -o ghs_built_v_3ss.csv
    """


if __name__ == "__main__":
    """
    # all - nres = res
    gdal_calc.py \
        -A GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0__THA.tif \
        -B GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0__THA.tif \
        --outfile="ghs_built_v_res_3ss__THA.tif" \
        --calc="A-B"

    cp GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0__THA.tif \
        ghs_built_v_nres_3ss__THA.tif
    """
    value_150ss_tif = "input/giri/THA/bem_5x5_valfis_res__THA.tif"
    volume_3ss_tif = "input/ghsl/THA/ghs_built_v_res_3ss__THA.tif"
    flood_1ss_tif = "input/footprints/JBA/Raster/TH_FLRF_ChaoPhraya2011_RD_01.tif"
    prefix = "res"
    main(value_150ss_tif, volume_3ss_tif, flood_1ss_tif, prefix)

    value_150ss_tif = "input/giri/THA/bem_5x5_valfis_nres__THA.tif"
    volume_3ss_tif = "input/ghsl/THA/ghs_built_v_nres_3ss__THA.tif"
    flood_1ss_tif = "input/footprints/JBA/Raster/TH_FLRF_ChaoPhraya2011_RD_01.tif"
    prefix = "nres"
    main(value_150ss_tif, volume_3ss_tif, flood_1ss_tif, prefix)
