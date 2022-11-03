import logging
from glob import glob

import geopandas as gpd
import pandas as pd
from tqdm import tqdm


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
logging.info("Start")

def get_damages():
    damage_fnames = glob("results/grid_damage_*")
    dfs = []
    for fname in tqdm(damage_fnames):
        df = pd.read_parquet(fname)
        dfs.append(df)
    logging.info("Done reading damages")

    damages = pd.concat(dfs)
    logging.info("Done concat")
    return damages

damages = get_damages().reset_index().groupby('edge_id').sum()

grid = gpd.read_parquet("results/grid.geoparquet").set_index('edge_id')
logging.info("Done reading grid")

grid_with_damages = grid.join(damages, validate='one_to_one')
logging.info("Done join")

grid_with_damages.to_parquet("results/grid_damages.geoparquet")
logging.info("Done parquet")

grid_with_damages[[c for c in grid_with_damages.columns if "ead" in c]+["source", "geometry"]].to_file("results/grid_damages.gpkg", driver="GPKG")
logging.info("Done GPKG")
