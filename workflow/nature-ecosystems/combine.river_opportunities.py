import sys

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.features

from open_gira.io import write_empty_frames


def query_ds_points(features, grid_points):
    """Runs in around 10-20s for 1000-100,000 features"""
    # check intersection
    joined = (
        features[["feature_id", "geometry"]]
        .sjoin_nearest(grid_points)
        .drop(columns=["geometry", "index_right"])
    )

    # mean values (generally per ha)
    return joined.groupby("feature_id").agg(np.nanmean)


if __name__ == "__main__":
    slice_n = snakemake.wildcards.SLICE_SLUG.replace("slice-", "")
    grid_points = gpd.read_parquet(
        snakemake.input.nbs_stack,
        columns=[
            "biodiversity_benefit",
            "carbon_benefit_t_per_ha",
            "planting_cost_usd_per_ha",
            "regen_cost_usd_per_ha",
            "geometry",
        ],
    )
    basins = gpd.read_parquet(
        snakemake.input.hydrobasins_adm,
        columns=["HYBAS_ID", "GID_0", "GID_1", "GID_2", "geometry"],
    )

    # Read opportunities
    with rasterio.open(snakemake.input.tree_potential_tif) as src:
        features = rasterio.features.dataset_features(src, bidx=1, geographic=False)
        options = gpd.GeoDataFrame.from_features(features)

    if options.empty:
        write_empty_frames(snakemake.output.parquet)
        sys.exit()

    options = options.drop(columns=["filename", "val"]).set_crs(src.crs)

    # Overlay Hydrobasins
    options = options.overlay(basins)

    # Set ID, area properties
    options["feature_id"] = [f"{slice_n}_{n}" for n in range(len(options))]
    options["area_m2"] = options.geometry.to_crs("ESRI:54009").area
    options["area_ha"] = options["area_m2"] * 1e-4

    # Read stacked points (costs and benefits)
    options_grid = query_ds_points(options, grid_points)
    options_attributed = options.set_index("feature_id").join(options_grid)
    options_attributed.to_parquet(snakemake.output.parquet)
