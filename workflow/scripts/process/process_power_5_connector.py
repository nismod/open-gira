"""Documents the connections between boxes"""


from process_power_functions import adj
from importing_modules import *

network_paths = glob.glob(
    os.path.join("data", "processed", "all_boxes", "box_*", "network_box_*.gpkg")
)  # finds the network_{box_id}.gpkg files that exist in all_boxes

for network_path in tqdm(
    network_paths, desc="connecting all boxes", total=len(network_paths)
):
    s = time.time()
    # extract box numbers
    filename = network_path[network_path.find("network_box_") :]
    box_id = filename[filename.find("box_") : filename.find(".gpkg")]
    id = int(box_id[4:])
    examine = adj(id)

    portal_dict = {}  # master dict
    gdf_idx = gpd.read_file(network_path, layer="nodes")
    if list(gdf_idx.id) != [None]:  # check not dummy

        for id_ex in examine:
            portal_lst_ex = []  # neighbour list
            path_test = os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"box_{id_ex}",
                f"network_box_{id_ex}.gpkg",
            )
            if os.path.exists(path_test):
                gdf_ex = gpd.read_file(path_test, layer="nodes")
                if list(gdf_ex.id) != [None]:  # check not dummy
                    # points_ex = gdf_ex.overlay(gdf_idx, how='intersection')  # find overlaps
                    # points_idx = gdf_idx.overlay(gdf_ex, how='intersection')  # find overlaps
                    duplicates = gdf_ex.merge(
                        gdf_idx,
                        left_on="geometry",
                        right_on="geometry",
                        suffixes=("_ex", "_idx"),
                    )
                    # update_dict = {i: j for i, j in zip(duplicates.id_idx, duplicates.id_ex)}  # add {examined_box_point:other_box_point}
                    update_lst = [
                        [
                            f"edge_X_{from_box}__{to_box}",
                            from_source_id,
                            f"cb_{from_box}__{to_box}",
                            "transmission",
                            from_idx,
                            to_ex,
                            str(LineString([fromto_geom, fromto_geom])),
                            f"{from_idx}__{to_ex}",
                        ]
                        for from_idx, to_ex, from_box, to_box, from_source_id, fromto_geom in zip(
                            duplicates.id_idx,
                            duplicates.id_ex,
                            duplicates.box_id_idx,
                            duplicates.box_id_ex,
                            duplicates.source_id_idx,
                            duplicates.geometry,
                        )
                    ]
                    # update_dict = {i: j for i, j in zip(duplicates.id_idx, duplicates.id_ex)}
                    # portal_lst_ex+=update_lst

                    portal_dict.update(
                        {f"box_{id_ex}": update_lst}
                    )  # update master dict with neighbour dict
    with open(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"connector_{box_id}.txt"
        ),
        "w",
    ) as file_ex:
        json.dump(portal_dict, file_ex)
