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
    slice_n = snakemake.wildcards.SLICE_SLUG.replace("slice-", "")  # noqa: F821
    grid_points = gpd.read_parquet(
        snakemake.input.nbs_stack,  # noqa: F821
        columns=[
            "biodiversity_benefit",
            "carbon_benefit_t_per_ha",
            "mangrove_planting_cost_usd_per_ha",
            "mangrove_regen_cost_usd_per_ha",
            "geometry",
        ],
    ).rename(
        columns={
            "mangrove_planting_cost_usd_per_ha": "planting_cost_usd_per_ha",
            "mangrove_regen_cost_usd_per_ha": "regen_cost_usd_per_ha",
        }
    )
    basins = gpd.read_parquet(
        snakemake.input.hydrobasins_adm,  # noqa: F821
        columns=["HYBAS_ID", "GID_0", "GID_1", "GID_2", "geometry"],
    ).reset_index()

    # Read opportunities
    with rasterio.open(snakemake.input.mangrove_potential_tif) as src:  # noqa: F821
        features = rasterio.features.dataset_features(src, bidx=1, geographic=False)
        options = gpd.GeoDataFrame.from_features(features)

    if options.empty:
        write_empty_frames(snakemake.output.parquet)  # noqa: F821
        sys.exit()

    options = options.query("val != 0").drop(columns="filename").set_crs(src.crs)
    labels = {
        1: "accreting",  # Accreting (expanding) shorelines. Ideal
        2: "retreating",  # Static to moderate retreating shorelines (sub ideal)
        3: "retreating_fast",  # Fast retreating shorelines.
        # ..might be worth considering if coastal managment is applied to block fill an otherwise appropriate area
    }
    options["option_shoreline"] = options.val.map(labels)
    options.drop(columns="val", inplace=True)

    # NOTE: for mangroves, join to nearest hydrobasin (don't split)
    options = options.sjoin_nearest(basins).drop(columns="index_right")

    # Set ID, area properties
    options["feature_id"] = [f"{slice_n}_{n}" for n in range(len(options))]
    options["area_m2"] = options.geometry.to_crs("ESRI:54009").area
    options["area_ha"] = options["area_m2"] * 1e-4

    # Read stacked points (costs and benefits)
    options_grid = query_ds_points(options, grid_points)
    options_attributed = options.set_index("feature_id").join(options_grid)
    options_attributed.to_parquet(snakemake.output.parquet)  # noqa: F821
