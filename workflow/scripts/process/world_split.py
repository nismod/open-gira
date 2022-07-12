"""Splits the world into boxes of length and height boxlen (input parameter)"""


# N.B. this is evil
# TODO: remove importing_modules and have script specific imports
from importing_modules import *

from collections import defaultdict

from shapely.strtree import STRtree


try:
    boxlen = snakemake.params["boxlen_value"]
    output_dir = snakemake.params["output_dir"]
except:
    output_dir = sys.argv[1]
    boxlen = sys.argv[2]


boxlen = float(boxlen)


if __name__ == "__main__":

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
    with fiona.open(
        os.path.join(output_dir, "input", "admin-boundaries", f"gadm36_levels.gpkg"),
        "r",
        layer=0,
    ) as src_code:
        code_geoms = []
        code_GIDs = []
        for feature in src_code:
            code_geoms.append(shape(feature["geometry"]))
            code_GIDs.append(feature["properties"]["GID_0"])
        print("create dataframe")
        countries = gpd.GeoDataFrame({"geometry": code_geoms, "code": code_GIDs})

    print(f"country length: {len(countries)}")
    print("Find countries intersecting with each grid cell...")
    # build a spatial index from the grid polygons
    tree = STRtree(grid.geometry)

    # when querying the STRtree, we are returned intersecting grid cell polygons
    # we really want their ids, so we need a mapping to recover those
    # shapely polygons are not hashable (at least not until shapely 2.0)
    # so use something that is hashable, their WKT representation (a string)
    box_wkt_to_id: dict[str, str] = {box.geometry.wkt: box.box_id for box in grid.itertuples()}

    # "box_n" -> ['GBR', 'FRA']
    countries_by_box: dict[str, list[str]] = defaultdict(list)

    for country in countries.itertuples():
        # query the tree for boxes which overlap the country
        intersecting_box_ids: list[str] = [box_wkt_to_id[box.wkt] for box in tree.query(country.geometry)]

        # add the country code to every boxes's list
        [countries_by_box[box_id].append(country.code) for box_id in intersecting_box_ids]

    # post-processing to match Max's previous output
    # 1) assign None to empty boxes
    for box in grid.itertuples():
        if box.box_id not in countries_by_box:
            countries_by_box[box.box_id] = None
    # 2) sort the elements by their box index
    countries_by_box = dict(sorted(countries_by_box.items(), key=lambda item: int(item[0].split("_")[-1])))

    print("writing metadata...")
    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.txt"), "w"
    ) as filejson:
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
        json.dump(info, filejson)

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
    grid.to_file(
        os.path.join(output_dir, "power_processed", "world_boxes.gpkg"), driver="GPKG"
    )
