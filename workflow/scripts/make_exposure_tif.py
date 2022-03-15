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
import os
import numpy as np
import sys

if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    try:
        db_file = snakemake.input['geoparquet']
        hazard = snakemake.input['hazard'][0]
        output_path = snakemake.output[0]
        threshold = snakemake.config['exposure_threshold']
    except NameError:
        print(sys.argv)
        (
            db_file,
            hazard,
            output_path,
            threshold
        ) = sys.argv[1:]
        # db = '../../results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet'
        # hazard = '../../results/input/hazard-aqueduct-river/tanzania-mini/inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif'
        # output_path = '../../results/test.tif'
        # threshold = 0.25
    output_path = os.path.join(output_path, f"exposure_{os.path.basename(hazard)}")
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
                output_path,
                'w',
                driver=h_file.driver,
                height=h_file.height,
                width=h_file.width,
                count=1,
                dtype=db['length_km'].dtype,
                crs=h_file.crs,
                transform=h_file.transform
        ) as output_file:
            # Zero the data and then fill in non-missing values
            # If a faster approach is needed we can probably find one using clever masking
            # tricks to write everything at once
            data = np.zeros((h_file.height, h_file.width), db['length_km'].dtype)
            for i, r in db.iterrows():
                data[r['cell_y'], r['cell_x']] = r['length_km']
            output_file.write(data, 1)
