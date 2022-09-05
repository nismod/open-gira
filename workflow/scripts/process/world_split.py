"""
Splits the world into boxes of length and height boxlen (input parameter)
"""


# N.B. this is evil
# TODO: remove importing_modules and have script specific imports
from importing_modules import *


if __name__ == "__main__":

    try:
        admin_data_path = snakemake.input["admin_data"]
        boxlen = snakemake.params["boxlen_value"]
        output_dir = snakemake.params["output_dir"]
        global_metadata_path = snakemake.output["global_metadata_path"]
        global_boxes_path = snakemake.output["global_boxes_path"]
    except:
        admin_data_path = sys.argv[1]
        output_dir = sys.argv[2]
        boxlen = sys.argv[3]
        global_metadata_path = sys.argv[4]
        global_boxes_path = sys.argv[5]

    boxlen = float(boxlen)

    assert 180 % boxlen == 0  # ensure divisibility

    lat_centroids = np.arange(90 - boxlen / 2, -90, -boxlen)
    lon_centroids = np.arange(-180 + boxlen / 2, 180, boxlen)

    points = []

    for lat in lat_centroids:  # left to right and down
        for lon in lon_centroids:
            points.append(Point(lon, lat))

    points_gdf = gpd.GeoDataFrame({"geometry": points})
    grid = points_gdf.buffer(boxlen / 2, cap_style=3)
    grid = gpd.GeoDataFrame(
        {
            "geometry": grid,
            "box_id": [f"box_{idx_}" for idx_ in range(len(grid))],
        },
        crs="EPSG:4326",
    )

    print("loading country layer=0")
    countries = gpd.read_file(admin_data_path, layer=0).drop(["NAME_0"], axis="columns")
    countries = countries.rename({"GID_0": "code"}, axis="columns")

    print(f"country length: {len(countries)}")
    print("Find countries intersecting with each grid cell...")

    # e.g. "box_284" -> ['GBR', 'FRA']
    countries_by_box: dict[str, list[str]] = {}
    for box_id, box_countries in gpd.sjoin(countries, grid).groupby('box_id'):
        countries_by_box[box_id] = list(box_countries.code)

    # post-processing to match Max's previous output
    # 1) assign None to empty boxes
    for box in grid.itertuples():
        if box.box_id not in countries_by_box:
            countries_by_box[box.box_id] = None
    # 2) sort the elements by their box index
    countries_by_box = dict(sorted(countries_by_box.items(), key=lambda item: int(item[0].split("_")[-1])))

    print("writing metadata...")
    with open(global_metadata_path, "w") as filejson:
        lon_min, lat_min, lon_max, lat_max = grid.bounds.values[0]
        info = {
            "boxlen": boxlen,
            "lon_min": lon_min,
            "lat_min": lat_min,
            "lon_max": lon_max,
            "lat_max": lat_max,
            "num_cols": len(lon_centroids),
            "num_rows": len(lat_centroids),
            "tot_boxes": int(len(lon_centroids) * len(lat_centroids)),
            "box_country_dict": countries_by_box,
        }
        json.dump(info, filejson, indent=2, sort_keys=True)

    print("creating box folders")
    for box_id in grid["box_id"]:
        all_boxes_path = os.path.join(
            output_dir, "power_processed", "all_boxes", f"{box_id}"
        )
        if not os.path.exists(all_boxes_path):
            os.makedirs(all_boxes_path)

    # for 5 degree boxes, this means writing ~2600 geopackage files
    # TODO: investigate if we really need to write these like this
    for box_id, box in tqdm(
        grid.groupby("box_id"), desc="saving each box", total=len(grid)
    ):  # separately so that world_boxes.gpkg can be opened on QGIS without having to rerun snakemake
        box.to_file(
            os.path.join(
                output_dir,
                "power_processed",
                "all_boxes",
                box_id,
                f"geom_{box_id}.gpkg",
            ),
            driver="GPKG",
        )

    print("writing full to gpkg")
    grid.to_file(global_boxes_path, driver="GPKG")
