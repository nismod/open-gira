"""
Download and merge JRC river flood return period maps

https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_rp50y-tif

The global river flood hazard maps are a gridded data set representing
inundation along the river network, for seven different flood return periods
(from 1-in-10-years to 1-in-500-years). The input river flow data for the new
maps are produced by means of the open-source hydrological model LISFLOOD, while
inundation simulations are performed with the hydrodynamic model LISFLOOD-FP.
The extent comprises the entire world with the exception of Greenland and
Antarctica and small islands with river basins smaller than 500km2.

Cell values indicate water depth (in m). The maps can be used to assess the
exposure of population and economic assets to river floods, and to perform flood
risk assessments. The dataset is created as part of the Copernicus Emergency
Management Service. NOTE: this dataset is not an official flood hazard map (for
details and limitations please refer to related publications).

Dataset is tiled. This script first downloads the tiles and then merges them
into a single global GeoTiff

Citation:

Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Dottori, Francesco; Neal,
Jeffrey; Prudhomme, Christel; Salamon, Peter (2024): Global river flood hazard
maps. European Commission, Joint Research Centre (JRC) [Dataset] PID:
http://data.europa.eu/89h/jrc-floods-floodmapgl_rp50y-tif
"""

import os
import requests
import logging

import geopandas as gpd
import rasterio
from osgeo import gdal

if __name__ == "__main__":
    try:
        raw_folder: str = snakemake.input["raw_folder"]
        RP10_file = snakemake.output["RP10"]
        RP20_file = snakemake.output["RP20"]
        RP50_file = snakemake.output["RP50"]
        RP75_file = snakemake.output["RP75"]
        RP100_file = snakemake.output["RP100"]
        RP200_file = snakemake.output["RP200"]
        RP500_file = snakemake.output["RP500"]
        # base_directory = snakemake.wildcard.OUTPUT_DIR
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
# What RPs are we working with?
RPS = ["10", "20", "50", "75", "100", "200", "500"]
RP_files = [RP10_file, RP20_file, RP50_file, RP75_file, RP100_file, RP200_file, RP500_file]

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Preparing files")
# os.makedirs(raw_folder, exist_ok=True)
# os.makedirs(merged_folder, exist_ok=True)
VRT_FILE = os.path.join(raw_folder, "temp_vrt.vrt") # temp virtual raster for merging

logging.info("Loading tile extents from GeoJSON")
BASE_URL = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/" # JRC flood data base URL
geojson_path = os.path.join(BASE_URL, "tile_extents.geojson")
tiles = gpd.read_file(geojson_path)

for idx, RP in enumerate(RPS):
    logging.info(f"Working on RP {RP}.")
    raster_files = [] # List to store file paths for merging
    logging.info("Getting tile info.")
    for _, row in tiles.iterrows():
        tile_name = row["name"]  # e.g., "N80_W170"
        tile_id = row["id"]
        file_name = f"ID{tile_id}_{tile_name}_RP{RP}_depth.tif"
        file_path = f"{raw_folder}/{file_name}"

        # Can download again here if for some reason file is not already downloaded
        if not os.path.exists(file_path):
            print(f"Downloading {file_name}...")
            file_url = f"{BASE_URL}RP{RP}/{file_name}"
            response = requests.get(file_url, stream=True)
            if response.status_code == 200:
                with open(file_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=1024):
                        f.write(chunk)
                print(f"Saved: {file_path}")
            else:
                print(f"Failed to download {file_name}")

        raster_files.append(file_path)

    logging.info("Merging files into one global raster.")
    if raster_files:
        logging.info("Creating VRT file...")
        gdal.BuildVRT(VRT_FILE, raster_files, options=gdal.BuildVRTOptions(resampleAlg='nearest'))
    
    logging.info('Converting VRT to a merged GeoTiff...')
    MERGED_OUTPUT = RP_files[idx] # read filepath at same index of RP we are working on
    gdal.Translate(MERGED_OUTPUT, VRT_FILE, options=gdal.TranslateOptions(format="GTiff", creationOptions=["COMPRESS=LZW", "BIGTIFF=YES"]))
    
    if os.path.exists(VRT_FILE):
        logging.info('Deleting VRT file...')
        os.remove(VRT_FILE)

logging.info("Done")