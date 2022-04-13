"""
Creates the network from the plants and targets data

If box has:
no targets, no plants, no edges -> dummy file for all
no targets, no plants, edges -> dummy file for targets and plants, edges file
no targets, plants, no edges -> dummy target file, plants file but connected=False, dummy edges file
no targets, plants, edges -> dummy target file, plants file, edges file
targets, no plants, no edges -> targets file but connected=False, dummy plants file, dummy edges file
targets, plants, no edges -> targets file but connected=False, plants file but connected=False, dummy edges file
targets, no plants, edges -> targets file, dummy plants file, edges file
targets, plants, edges -> all files
"""

from importing_modules import *

try:
    box_id = snakemake.params["box_id"]
    output_dir = snakemake.params['output_dir']
except:
    output_dir = sys.argv[1]
    box_id = sys.argv[2]


def timer(s):
    print("timer : ", round((time.time() - s) / 60, 2), " mins\n")


#%%
if __name__ == "__main__":
    edges_empty = False
    plants_empty = False
    targets_empty = False

    print(f"{box_id}: opening files")
    start = time.time()
    plantsfile = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id, f"powerplants_{box_id}.csv"
    )
    targetsfile = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id, f"targets_{box_id}.csv"
    )

    plants = pd.read_csv(plantsfile)
    targets = pd.read_csv(targetsfile)

    plants["geometry"] = plants["geometry"].apply(wkt.loads)
    plants = gpd.GeoDataFrame(plants)

    targets["geometry"] = targets["geometry"].apply(wkt.loads)
    targets = gpd.GeoDataFrame(targets)

    plants = plants.reset_index().copy()
    plants["id"] = [f"source_{i}_{box_id}" for i in range(len(plants))]
    plant_cols = ["id", "source_id", "capacity_mw", "type", "box_id", "geometry"]
    plants = plants[plant_cols]
    plants[
        "connected"
    ] = False  # bool, will be set to true if edges available to be connected to
    print(f"{box_id}: writing to plants_{box_id}.gpkg")
    plants_file_name = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id, f"plants_{box_id}.gpkg"
    )
    if len(plants) == 0:  # if empty, write dummy
        plants_empty_df = gpd.GeoDataFrame(columns=plants.columns)
        plants_empty_df.loc[0, :] = [None] * len(plants.columns)
        plants_empty_df["box_id"] = box_id
        plants_empty_df.to_file(plants_file_name, driver="GPKG")
        plants_empty = True
    else:
        plants.to_file(plants_file_name, driver="GPKG")

    targets = targets.reset_index().copy()
    targets["id"] = [f"target_{i}_{box_id}" for i in range(len(targets))]

    # timer(start)

    print(f"{box_id}: writing to targets_{box_id}.gpkg")

    targets_file_name = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id, f"targets_{box_id}.gpkg"
    )
    target_cols = list(targets.columns)
    target_cols.remove("centroid")
    if len(targets) == 0:  # if empty, write dummy file
        targets_empty_df = gpd.GeoDataFrame(columns=target_cols)
        targets_empty_df.loc[0, :] = [None] * len(target_cols)
        targets_empty_df["box_id"] = box_id
        targets_empty_df.to_file(targets_file_name, driver="GPKG")
        targets_empty = True
    else:
        targets.drop(columns=["centroid"]).to_file(targets_file_name, driver="GPKG")
    # timer(start)

    # print("processing targets GeoData")
    targets = targets[["id", "type", "box_id", "centroid"]].rename(
        columns={"centroid": "geometry"}
    )
    # timer(start)

    # Combine to nodes
    # print("combining sources and sinks")
    nodes = plants.append(targets).reset_index(drop=True)  # sources and targets located
    # timer(start)

    # Edges
    # print("getting gridfinder lines")

    edges = gpd.read_file(
        os.path.join(
            output_dir, "power_processed", "all_boxes", f"{box_id}", f"gridfinder_{box_id}.gpkg"
        )
    )

    if len(edges) == 1 and edges["geometry"].any() == None:  # catch empty
        edges = gpd.GeoDataFrame(columns=edges.columns)
        edges_empty = True

    # timer(start)

    # print("processing lines")
    edges["type"] = "transmission"
    edges["id"] = edges.reset_index()["index"].apply(
        lambda i: f"edge_{i}_{box_id}"
    )  # changed to edges, was targets (?)

    # timer(start)

    # print("processing edges GeoData")
    edges = edges[["id", "source_id", "box_id", "type", "geometry"]]
    # timer(start)

    # Process network
    print(f"{box_id}: creating network")
    network = snkit.network.Network(nodes, edges)
    # timer(start)

    # fix str when should be Point(# #)
    network.nodes["geometry"] = [
        sw.loads(x) if type(x) == str else x for x in network.nodes["geometry"]
    ]
    network.edges["geometry"] = [
        sw.loads(x) if type(x) == str else x for x in network.edges["geometry"]
    ]

    network.nodes["connected"] = False  # default

    if not edges_empty:

        # Connect power plants
        # print("connecting powerplants")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
            warnings.simplefilter(action="ignore", category=FutureWarning)

            # print("processing network")
            network = snkit.network.split_multilinestrings(network)
            # timer(start)

            geod = Geod(ellps="WGS84")
            edge_limit = 20_000  # meters   # TODO 200_000 vs 20_000

            network = snkit.network.link_nodes_to_nearest_edge(
                network,
                lambda node, edge: node.type == "source"
                and geod.geometry_length(edge.geometry) < edge_limit,
            )

        network.nodes.loc[network.nodes.id.isnull(), "type"] = "conn_source"
        # timer(start)

        # print("resetting indices")
        network.nodes["id"] = network.nodes.reset_index().apply(
            lambda row: f"conn_source_{row['index']}_{box_id}"
            if type(row.id) is float
            else row.id,
            axis=1,
        )

        # Connect targets
        # print("connecting targets")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
            warnings.simplefilter(action="ignore", category=FutureWarning)

            network = snkit.network.link_nodes_to_nearest_edge(
                network,
                lambda node, edge: node.type == "target"
                and geod.geometry_length(edge.geometry) < edge_limit,
            )
        network.nodes.loc[network.nodes.id.isnull(), "type"] = "conn_target"
        # timer(start)

        # print("resetting indices")
        network.nodes["id"] = network.nodes.reset_index().apply(
            lambda row: f"conn_target_{row['index']}_{box_id}"
            if type(row.id) is float
            else row.id,
            axis=1,
        )
        # timer(start)

        # Add nodes at line endpoints
        # print("adding nodes at line endpoints")
        network = snkit.network.add_endpoints(network)
        # timer(start)

        # print("processing network")
        network.nodes.loc[network.nodes.id.isnull(), "type"] = "intermediate"
        # timer(start)

        # print("resetting indices")
        network.nodes["id"] = network.nodes.reset_index().apply(
            lambda row: f"intermediate_{row['index']}_{box_id}"
            if type(row.id) is float
            else row.id,
            axis=1,
        )
        # timer(start)

        # add from/to ids
        # print("adding from/to ids")
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            network = snkit.network.add_topology(network, id_col="id")
        # timer(start)

        network.edges["box_id"] = box_id

        network.nodes["box_id"] = box_id
        network.nodes[
            "connected"
        ] = True  # since can be connected to network (edges available)
    else:  # empty edges
        network.edges = gpd.GeoDataFrame(columns=network.edges.columns)
        network.edges.loc[0, :] = [None] * len(network.edges.columns)
        network.edges["box_id"] = box_id

    # output
    out_fname = os.path.join(
        output_dir, "power_processed", "all_boxes", box_id, f"network_{box_id}.gpkg"
    )
    print("writing edges to ", out_fname)
    network.edges.to_file(out_fname, layer="edges", driver="GPKG")
    # timer(start)

    if targets_empty and plants_empty:
        network.nodes = gpd.GeoDataFrame(columns=network.nodes.columns)
        network.nodes.loc[0, :] = [None] * len(network.nodes.columns)
        network.nodes["box_id"] = box_id

    print("writing nodes to ", out_fname)
    network.nodes.to_file(out_fname, layer="nodes", driver="GPKG")
    # timer(start)

    print(
        "\nnote that ShapelyDepreciationWarnings (shapely) and FutureWarning (geopandas) were supressed. Newer pandas versions could cause issues (1.1.5 works)."
    )  # remove if fixed
    end = time.time()
    print(f"\n{box_id}: Time to run file: {(end - start)/60:0.2f} minutes")
