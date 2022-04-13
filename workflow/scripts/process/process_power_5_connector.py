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
    box_id = os.path.basename(network_path)[8:-5]
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
                    duplicates = gdf_ex.merge(
                        gdf_idx,
                        left_on="geometry",
                        right_on="geometry",
                        suffixes=("_ex", "_idx"),
                    )
                    update_lst = [
                        {
                            "edge_id": f"edge_X_{from_box}__{to_box}",
                            "source_id": from_source_id,
                            "link": f"cb_{from_box}__{to_box}",
                            "type": "transmission",
                            "from_id": from_idx,  # connection in box_idx
                            "to_id": to_ex,  # connection in box next to box_idx
                            "geometry": str(LineString([fromto_geom, fromto_geom])),
                        }
                        for from_idx, to_ex, from_box, to_box, from_source_id, fromto_geom in zip(
                            duplicates.id_idx,
                            duplicates.id_ex,
                            duplicates.box_id_idx,
                            duplicates.box_id_ex,
                            duplicates.source_id_idx,
                            duplicates.geometry,
                        )
                    ]

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
