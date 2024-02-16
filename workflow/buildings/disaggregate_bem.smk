"""
Upscale (disaggregate) building exposure layers according to built volume
"""

rule disaggregate_bem:
    input:
        bem_res=rules.download_giri_bem.output.res,
        ghsl_res="{OUTPUT_DIR}/input/ghsl/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        bem_nres=rules.download_giri_bem.output.nres,
        ghsl_nres="{OUTPUT_DIR}/input/ghsl/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
    output:
        bem_res_3ss="{OUTPUT_DIR}/buildings/building_exposure_res_3ss.tif",
        bem_nres_3ss="{OUTPUT_DIR}/buildings/building_exposure_nres_3ss.tif",
    shell:
        import rasterio
        from rasterio.warp import reproject, Resampling

        def disaggregate(value_150ss_ds, volume_3ss_ds):
            # BEM gives ~5x5km value in USD => value_150ss
            value_150ss = value_150ss_ds.read(1)
            print("Read value", value_150ss.shape, value_150ss.dtype)

            # GHSL gives ~100m  volume in m3 => volume_3ss
            volume_3ss = volume_3ss_ds.read(1)
            print("Read volume", volume_3ss.shape, volume_3ss.dtype)


            # Resample GHSL to ~5x5km, "sum" => volume_150ss
            factor = 0.02
            volume_150ss = volume_3ss_ds.read(
                out_shape=(
                    volume_3ss.count,
                    int(volume_3ss.height * factor),
                    int(volume_3ss.width * factor)
                ),
                resampling=Resampling.sum
            )
            print("Read volume coarse", volume_150ss.shape, volume_150ss.dtype)


            # Calculate (coarse) value per volume
            value_per_volume_150ss = value_150ss / volume_150ss

            # Resample to fine-scale value per volume, "nearest"
            with rasterio.Env():
                value_per_volume_3ss = np.zeros(volume_3ss.shape, np.float64)
                reproject(
                    value_per_volume_150ss,
                    value_per_volume_3ss,
                    src_transform=volume_150ss.transform,
                    src_crs=volume_150ss.crs,
                    dst_transform=volume_3ss_ds.transform,
                    dst_crs=volume_3ss_ds.crs,
                    resampling=Resampling.nearest)

            # Calculate fine-scale value
            value_3ss = value_per_volume_3ss * volume_3ss

            return value_3ss

        def write_raster_like(data, filename, template_ds):
            with rasterio.open(
                    filename,
                    'w',
                    driver='GTiff',
                    width=template_ds.width,
                    height=template_ds.height,
                    count=1,
                    dtype=np.float64,
                    nodata=0,
                    transform=template_ds.transform,
                    crs=template_ds.transform) as output_ds:
                output_ds.write(data, indexes=1)

        # Residential
        res_value_150ss_ds = rasterio.open(input.bem_res)
        res_volume_3ss_ds = rasterio.open(input.ghsl_res)
        res_value_3ss = disaggregate(res_value_150ss_ds, res_volume_3ss_ds)
        write_raster_like(res_value_3ss, output.bem_res_3ss, res_volume_3ss_ds)
        res_value_150ss_ds.close()
        res_volume_3ss_ds.close()

        # Non-residential
        nres_value_150ss_ds = rasterio.open(input.bem_res)
        nres_volume_3ss_ds = rasterio.open(input.ghsl_res)
        nres_value_3ss = disaggregate(nres_value_150ss_ds, nres_volume_3ss_ds)
        write_raster_like(nres_value_3ss, output.bem_nres_3ss, nres_volume_3ss_ds)
        nres_value_150ss_ds.close()
        nres_volume_3ss_ds.close()



        # Then to assess flood damage

        # JBA footprint gives ~30m depth in m => depth_30m
        # Resample value_3ss to ~30m,  "nearest", divide by 9 => value_30m
        # Apply damage curve to depth_30m => damage_fraction_30m
        # Calculate (damage_fraction_30m * value_30m) => damage_30m
