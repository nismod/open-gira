"""Splits the world into boxes of length and height boxlen (input parameter)"""


from importing_modules import *

try:
    boxlen = snakemake.params["boxlen_value"]
    output_dir = snakemake.params['output_dir']
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
    gdf_area = points_gdf.buffer(boxlen / 2, cap_style=3)
    gdf_area = gpd.GeoDataFrame(
        {
            "geometry": gdf_area,
            "box_id": [f"box_{idx_}" for idx_ in range(len(gdf_area))],
        },
        crs="EPSG:4326",
    )

    print("loading country layer=0")
    with fiona.open(
        os.path.join(output_dir, "input", "adminboundaries", f"gadm36_levels.gpkg"), "r", layer=0
    ) as src_code:
        code_geoms = []
        code_GIDs = []
        for feature in src_code:
            code_geoms.append(shape(feature["geometry"]))
            code_GIDs.append(feature["properties"]["GID_0"])
        print("create dataframe")
        code_geoms_gpd = gpd.GeoDataFrame({"geometry": code_geoms, "code": code_GIDs})

    print(f"country length: {len(code_geoms_gpd)}")
    print("Country for each box...")
    box_country_dict = {}
    for jj, box in tqdm(
        enumerate(gdf_area["geometry"]),
        desc="intersecting countries",
        total=len(gdf_area),
    ):  # find countries that are in the boxes
        countries_overlap = []
        for ii, country in enumerate(code_geoms_gpd["geometry"]):
            if box.intersects(country) == True:
                countries_overlap.append(code_geoms_gpd["code"].iloc[ii])
        if len(countries_overlap) == 0:
            countries_overlap = None
        box_country_dict[gdf_area["box_id"].iloc[jj]] = countries_overlap

    print("writing metadata...")
    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.txt"), "w"
    ) as filejson:
        lon_min, lat_min, lon_max, lat_max = gdf_area.bounds.values[0]
        info = {
            "boxlen": boxlen,
            "lon_min": lon_min,
            "lat_min": lat_min,
            "lon_max": lon_max,
            "lat_max": lat_max,
            "num_cols": len(lon_centroids),
            "num_rows": len(lat_centroids),
            "tot_boxes": int(len(lon_centroids) * len(lat_centroids)),
            "box_country_dict": box_country_dict,
        }
        json.dump(info, filejson)

    print("creating box folders")
    for box_id in gdf_area["box_id"]:
        all_boxes_path = os.path.join(output_dir, "power_processed", "all_boxes", f"{box_id}")
        if not os.path.exists(all_boxes_path):
            os.makedirs(all_boxes_path)

    for box_id, gdf_box in tqdm(
        gdf_area.groupby("box_id"), desc="saving each box", total=len(gdf_area)
    ):  # separately so that world_boxes.gpkg can be opened on QGIS without having to rerun snakemake
        gdf_box.to_file(
            os.path.join(
                output_dir, "power_processed", "all_boxes", box_id, f"geom_{box_id}.gpkg"
            ),
            driver="GPKG",
        )

    print("writing full to gpkg")
    gdf_area.to_file(
        os.path.join(output_dir, "power_processed", "world_boxes.gpkg"), driver="GPKG"
    )
