"""Creates unit grid at resolution of the return period maps."""


import numpy as np
import pandas as pd
import netCDF4 as nc4
import sys
import geopandas as gpd
import os
import time
import json
from shapely.geometry import box

from pathos.multiprocessing import ProcessPool, cpu_count

from tqdm import tqdm

region = sys.argv[1]

squarehalfwidth = (
    0.05  # this is half the width of the smallest unit in the return period maps
)


def t(num, t):
    print(str(num), time.time() - t)


def make_grid_points_nc2(box_id, region, ps):
    """Updated grid point maker, uses return period maps
    Performs manual overlay, remember to not do this later then."""

    fn = os.path.join(
        "data", "stormtracks", "fixed", f"STORM_FIXED_RETURN_PERIODS_{region}.nc"
    )
    ds = nc4.Dataset(fn)
    lons = np.array(ds["lon"])
    lats = np.array(ds["lat"])

    assert (
        abs(np.average(lons[1:] - lons[:-1]) - 2 * squarehalfwidth) < 1e-12
    )  # checks that the defined squarehalfwidth is correct
    assert (
        abs(np.average(lats[1:] - lats[:-1]) - 2 * squarehalfwidth) < 1e-12
    )  # see above

    box_gs = gpd.read_file(
        os.path.join("data", "processed", "all_boxes", box_id, f"geom_{box_id}.gpkg")
    )
    lon_min, lat_min, lon_max, lat_max = box_gs.bounds.values[0]

    lons = lons[lons > lon_min]
    lons = lons[lons < lon_max]
    lats = lats[lats > lat_min]
    lats = lats[lats < lat_max]

    assert len(lons) != 0
    assert len(lats) != 0

    point_df = pd.DataFrame()
    tot = len(lats) * len(lons)

    box_infrastructure = gpd.read_file(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"gridfinder_{box_id}.gpkg"
        )
    )
    containing_box_dict = {}

    unit = 0
    for lat in tqdm(lats, desc=ps + " Manual Overlay", total=len(lats)):
        # for lat in lats:
        for (
            lon
        ) in (
            lons
        ):  # another option is to save which geoms are within then use this later rather than overlay again

            box_infrastructure_indiv = box_infrastructure[
                box_infrastructure["geometry"].intersects(
                    box(
                        lon - squarehalfwidth,
                        lat - squarehalfwidth,
                        lon + squarehalfwidth,
                        lat + squarehalfwidth,
                    )
                )
            ]  # DataFrame of which gridfinder lines are within the unit

            if (
                len(box_infrastructure_indiv) != 0
            ):  # if there are elemetns of infrastructure within/intersects unit, keep, check this length, if 0 then nothing there, don't include, else include
                ID_point = f"{region}_{box_id}_{unit}"
                containing_box_dict[ID_point] = list(
                    box_infrastructure_indiv["source_id"]
                )  # add the elements which are in the unit
                point_df_add = pd.DataFrame(
                    {"longitude": [lon], "latitude": [lat], "ID_point": [ID_point]}
                )
                point_df = pd.concat([point_df, point_df_add])
            unit += 1
    try:
        gdf = gpd.GeoDataFrame(
            point_df, geometry=gpd.points_from_xy(point_df.longitude, point_df.latitude)
        )
        return gdf, containing_box_dict
    except:  # no points at all
        return [], containing_box_dict


def create_grid_box(box_id, idx, totboxes):

    ps = f"{idx}/{totboxes}"  # print statement
    print(ps)
    ### create grid
    grid_box_indiv, containing_box_dict = make_grid_points_nc2(box_id, region, ps)
    if len(grid_box_indiv) != 0:  # exclude empty boxes

        grid_box_indiv = grid_box_indiv.reset_index(drop=True)

        ### extract points inside
        grid_box_indiv["box_id"] = box_id

        grid_box_indiv["region"] = region

        return grid_box_indiv, containing_box_dict


if __name__ == "__main__":
    nodesuse = max(1, cpu_count() - 2)
    if "linux" not in sys.platform:
        nodesuse = 8

    with open(
        os.path.join("data", "intersection", "regions", f"{region}_boxes.txt"), "r"
    ) as src:
        box_ids = json.load(src)
    totboxes = [len(box_ids)] * len(box_ids)
    idx_bxs = list(range(len(totboxes)))

    pool_grid = ProcessPool(nodes=nodesuse)
    output = pool_grid.map(create_grid_box, box_ids, idx_bxs, totboxes)

    output_grid = [item[0] for item in output]  # extract dataframes
    grid_boxes = pd.concat(output_grid).reset_index(drop=True)

    grid_boxes_area_series = grid_boxes.geometry.buffer(squarehalfwidth, cap_style=3)
    grid_boxes.rename(columns={"geometry": "centroid"}, inplace=True)

    grid_boxes_area = grid_boxes.assign(geometry=grid_boxes_area_series)

    grid_boxes_area["centroid"] = grid_boxes_area["centroid"].astype(str)
    grid_boxes_area.to_file(
        os.path.join("data", "intersection", "regions", f"{region}_unit.gpkg"),
        driver="GPKG",
    )

    output_contains = [item[1] for item in output]  # extract dictionaries
    unit_contains = {}
    for dict_indiv in output_contains:  # join the dictionaries
        unit_contains.update(dict_indiv)
    with open(
        os.path.join("data", "intersection", "regions", f"{region}_unit_contains.txt"),
        "w",
    ) as writefile:
        json.dump(unit_contains, writefile)
