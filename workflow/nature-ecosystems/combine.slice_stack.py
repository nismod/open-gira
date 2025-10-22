from pathlib import Path

import pandas as pd
import geopandas as gpd
import xarray as xr


def stack(layers: pd.DataFrame, target_path: Path):
    layer_paths = layers.path.tolist()
    var = xr.Variable("key", layers.key.tolist())
    ds = (
        xr.concat(
            [
                xr.open_dataset(layer_path, engine="rasterio")
                for layer_path in layer_paths
            ],
            dim=var,
        )
        .squeeze("band", drop=True)
        .drop_vars("spatial_ref")
    )
    # Trade-off in chunk size vs number of files
    # (smaller chunks -> more files -> slower to write, unknown effect on reads)
    dsc = ds.chunk({"x": 10000, "y": 10000, "key": 1000})
    dsc.to_zarr(target_path, mode="w-")


def stack_to_points(grid_fname: Path, crs):
    grid_data = (
        xr.open_zarr(grid_fname)
        .to_dataframe()
        .reset_index()
        .pivot(index=["y", "x"], columns="key", values="band_data")
        .reset_index()
    )
    # as points
    grid_points = gpd.GeoDataFrame(
        data=grid_data.drop(columns=["y", "x"]),
        geometry=gpd.points_from_xy(grid_data.x, grid_data.y),
        crs=crs,
    )
    grid_points.to_parquet(grid_fname.replace("zarr", "parquet"))


if __name__ == "__main__":
    rasters = pd.DataFrame(
        {
            "key": [
                "biodiversity_benefit",
                "carbon_benefit_t_per_ha",
                "planting_cost_usd_per_ha",
                "regen_cost_usd_per_ha",
                "mangrove_planting_cost_usd_per_ha",
                "mangrove_regen_cost_usd_per_ha",
            ],
            "path": [
                snakemake.input.biodiversity_benefit_tif,  # noqa: F821
                snakemake.input.carbon_benefit_tif,  # noqa: F821
                snakemake.input.planting_cost_tif,  # noqa: F821
                snakemake.input.regeneration_cost_tif,  # noqa: F821
                snakemake.input.mangrove_planting_cost_tif,  # noqa: F821
                snakemake.input.mangrove_regeneration_cost_tif,  # noqa: F821
            ],
        }
    )
    stack(rasters, snakemake.output.zarr)  # noqa: F821
    stack_to_points(snakemake.output.zarr, "EPSG:4326")  # noqa: F821
