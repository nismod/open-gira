#!/usr/bin/env python
# coding: utf-8
"""
Create a raster file (.tif) from a .geoparquet file
by thresholding cells by flood depth and grading by
length of road within the cell.
"""
import logging
import geopandas as gp
import rasterio
from rasterio.enums import Resampling
import os
import numpy as np


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    # move to settings later
    epsg = 3857
    try:
        db_file = snakemake.input['geoparquet']
        hazard = snakemake.input['hazard'][0]
        coastline = snakemake.input['coastline']
        output_path = snakemake.output[0]
        opts_dict = snakemake.config['exposure_tifs']
    except NameError:
        # print(sys.argv)
        # (
        #     db_file,
        #     hazard,
        #     coastline,
        #     output_path,
        #     opts_dict
        # ) = sys.argv[1:]
        db_file = '../../results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet'
        hazard = '../../results/input/hazard-aqueduct-river/tanzania-mini/inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif'
        coastline = '../../results/input/coastlines/ne_10m_ocean/ne_10m_ocean.shp'
        output_path = '../../results/exposure/tanzania-mini/hazard-aqueduct-river'
        opts_dict = {}

    # Load up the options from the opts_dict
    threshold = opts_dict['exposure_threshold'] if 'exposure_threshold' in opts_dict.keys() else 0.5
    scaling_factor = opts_dict['scaling_factor'] if 'scaling_factor' in opts_dict.keys() else 0.1
    resampling_mode = opts_dict['resampling_mode'] if 'resampling_mode' in opts_dict.keys() else 'bilinear'
    resampling_mode = Resampling[resampling_mode]

    output_path = os.path.join(output_path, f"exposure_{os.path.basename(hazard)}")
    tmp_file_path = output_path.replace(".tif", ".tmp.tif")
    logging.info(f"Generating raster {output_path} from {db_file}")
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
        logging.info(f"Created directory {os.path.dirname(output_path)}")
    db = gp.read_parquet(db_file)
    col, ext = os.path.splitext(os.path.basename(hazard))
    db = (
        db.rename({col: 'hazard'}, axis=1)
            .query(f'hazard >= {float(threshold)}')
            .assign(
            cell_str=lambda x: x['cell_index'],
            cell_x=lambda x: [y[0] for y in x['cell_index']],
            cell_y=lambda x: [y[1] for y in x['cell_index']]
        )
            .astype({'cell_str': 'str'})  # wasn't working well without cell_str to group by
            .filter(['cell_x', 'cell_y', 'cell_str', 'length_km'])
            .groupby(['cell_str', 'cell_x', 'cell_y'])
            .agg({'length_km': 'sum'})
            .reset_index()  # ungroup
    )
    # Copy hazard file to the output destination, editing the band information
    with rasterio.open(hazard) as h_file:
        with rasterio.open(
                tmp_file_path,
                'w',
                driver=h_file.driver,
                height=h_file.height,
                width=h_file.width,
                count=1,
                dtype=db['length_km'].dtype,
                crs=h_file.crs,
                transform=h_file.transform
        ) as tmp_file:
            # Zero the data and then fill in non-missing values
            # If a faster approach is needed we can probably find one using clever masking
            # tricks to write everything at once
            data = np.zeros((h_file.height, h_file.width), db['length_km'].dtype)
            for i, r in db.iterrows():
                data[r['cell_y'], r['cell_x']] = r['length_km']
            tmp_file.write(data, 1)

        with rasterio.open(tmp_file_path) as tmp_file:
            scale_width = int(tmp_file.width * scaling_factor)
            scale_height = int(tmp_file.height * scaling_factor)

            # scale image transform
            transform = tmp_file.transform * tmp_file.transform.scale(
                (tmp_file.width / scale_width),
                (tmp_file.height / scale_height)
            )
            # Resampling
            with rasterio.open(
                    f"{output_path}",
                    'w',
                    driver=h_file.driver,
                    height=scale_height,
                    width=scale_width,
                    count=1,
                    dtype=db['length_km'].dtype,
                    crs=h_file.crs,
                    transform=transform
            ) as out_file:
                # resample data to target shape
                data = tmp_file.read(
                    out_shape=(
                        tmp_file.count,
                        scale_height,
                        scale_width
                    ),
                    resampling=resampling_mode
                )
                out_file.write(data)

    os.remove(tmp_file_path)
