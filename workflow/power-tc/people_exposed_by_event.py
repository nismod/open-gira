"""Calculate target population exposed to cyclone winds by event and threshold."""

import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

EVENT_CHUNK_SIZE = 128


def main() -> None:
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )

    logging.info("Reading target nodes")
    nodes: gpd.GeoDataFrame = gpd.read_parquet(snakemake.input.nodes)  # noqa: F821
    required_node_columns = {"asset_type", "population", "geometry"}
    missing_columns = required_node_columns - set(nodes.columns)
    if missing_columns:
        logging.warning(
            "Node file has no usable target schema; missing columns: %s. "
            "Writing an empty exposure table with the expected columns.",
            sorted(missing_columns),
        )
        targets = nodes.iloc[:0]
    else:
        targets = nodes.loc[nodes.asset_type == "target"].copy()

    if not targets.empty:
        targets = targets.to_crs("EPSG:4326")
        target_x = targets.geometry.x.to_numpy(dtype=float)
        target_y = targets.geometry.y.to_numpy(dtype=float)
        population = (
            pd.to_numeric(targets.population, errors="coerce")
            .fillna(0)
            .to_numpy(dtype=float)
        )
    else:
        target_x = target_y = population = np.array([], dtype=float)

    logging.info("Reading maximum wind fields")
    ds: xr.Dataset = xr.open_dataset(snakemake.input.wind_field)  # noqa: F821

    thresholds = list(snakemake.params.thresholds)  # noqa: F821
    threshold_columns = [f"{threshold:.1f}" for threshold in thresholds]
    event_ids = ds.event_id.values
    output_columns = ["event_id", *threshold_columns]

    if target_x.size and len(event_ids):
        # Multiple targets can share a raster cell. Sum their population once,
        # then sample and process each occupied cell only once per event.
        latitude_indices = ds.indexes["latitude"].get_indexer(target_y, method="nearest")
        longitude_indices = ds.indexes["longitude"].get_indexer(target_x, method="nearest")
        valid_targets = (latitude_indices >= 0) & (longitude_indices >= 0)
        cell_indices, target_to_cell = np.unique(
            np.column_stack(
                (latitude_indices[valid_targets], longitude_indices[valid_targets])
            ),
            axis=0,
            return_inverse=True,
        )
        cell_population = np.bincount(
            target_to_cell,
            weights=population[valid_targets],
            minlength=len(cell_indices),
        )
        cell_latitudes = ds.latitude.values[cell_indices[:, 0]]
        cell_longitudes = ds.longitude.values[cell_indices[:, 1]]

        logging.info(
            "Aggregated %d targets into %d occupied wind-field cells",
            int(valid_targets.sum()),
            len(cell_indices),
        )
        event_chunks = []

        logging.info("Calculating exposed population for %d events", len(event_ids))
        event_starts = range(0, len(event_ids), EVENT_CHUNK_SIZE)
        for start in tqdm(
            event_starts,
            total=(len(event_ids) + EVENT_CHUNK_SIZE - 1) // EVENT_CHUNK_SIZE,
            desc="Sampling wind by event batch",
            unit="batch",
        ):
            stop = min(start + EVENT_CHUNK_SIZE, len(event_ids))
            # NetCDF fancy indexing at many scattered points is very slow.
            # Read a contiguous event block, then gather occupied cells in memory.
            wind_grid_chunk = ds.max_wind_speed.isel(
                event_id=slice(start, stop)
            ).load()
            wind_at_cells = (
                wind_grid_chunk.sel(
                    latitude=xr.DataArray(cell_latitudes, dims="cell"),
                    longitude=xr.DataArray(cell_longitudes, dims="cell"),
                )
                .transpose("event_id", "cell")
                .to_numpy()
            )
            del wind_grid_chunk
            chunk = pd.DataFrame({"event_id": event_ids[start:stop].astype(str)})
            valid_wind = np.isfinite(wind_at_cells)
            for threshold, column in zip(thresholds, threshold_columns):
                exposed = valid_wind & (wind_at_cells >= threshold)
                chunk[column] = np.sum(
                    np.where(exposed, cell_population[np.newaxis, :], 0.0), axis=1
                )
            event_chunks.append(chunk)

        result = pd.concat(event_chunks, ignore_index=True)
    else:
        result = pd.DataFrame(columns=output_columns)
    output_path = snakemake.output.by_event  # noqa: F821
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    result.to_parquet(output_path, index=False)
    ds.close()


if __name__ == "__main__":
    main()
