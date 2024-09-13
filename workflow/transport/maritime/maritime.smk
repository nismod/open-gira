rule create_maritime_network:
    input:
        nodes = "{OUTPUT_DIR}/input/networks/maritime/nodes.gpq",
        edges_no_geom = "{OUTPUT_DIR}/input/networks/maritime/edges_by_cargo/maritime_base_network_general_cargo.pq",
        edges_visualisation = "{OUTPUT_DIR}/input/networks/maritime/edges.gpq",
    output:
        nodes = "{OUTPUT_DIR}/maritime_network/nodes.gpq",
        edges = "{OUTPUT_DIR}/maritime_network/edges.gpq",
    run:
        import igraph as ig
        import geopandas as gpd
        from shapely.geometry import Point
        from shapely.ops import linemerge
        from tqdm import tqdm

        from open_gira.network_creation import preprocess_maritime_network

        # possible cargo types = ("container", "dry_bulk", "general_cargo",  "roro", "tanker")
        # for now, just use 'general_cargo'
        maritime_nodes, maritime_edges_no_geom = preprocess_maritime_network(
            input.nodes,
            input.edges_no_geom
        )

        if config["study_country_iso_a3"] == "THA":
            # put Bangkok port in the right place...
            maritime_nodes.loc[maritime_nodes.name == "Bangkok_Thailand", "geometry"] = Point((100.5753, 13.7037))

        # %%
        # Jasper's maritime edges in 'edges_by_cargo' do not contain geometry
        # this is because the AIS data that they were derived from only contain origin and destination port, not route
        # this is a pain for visualisation, so we will create a geometry for each from `maritime_vis_edges`

        maritime_vis_edges = gpd.read_parquet(input.edges_visualisation)
        vis_graph = ig.Graph.DataFrame(maritime_vis_edges, directed=True, use_vids=False)

        maritime_edges = maritime_edges_no_geom.copy()
        change_of_port_mask = maritime_edges_no_geom.from_port != maritime_edges_no_geom.to_port
        port_pairs_to_generate_geom_for = maritime_edges_no_geom[change_of_port_mask]
        for index, row in tqdm(port_pairs_to_generate_geom_for.iterrows(), total=len(port_pairs_to_generate_geom_for)):
            edge_list = vis_graph.get_shortest_path(row.from_port, row.to_port, weights="distance", output="epath") 
            route_edges = maritime_vis_edges.iloc[edge_list]
            route_linestring = linemerge(list(route_edges.geometry))
            maritime_edges.loc[index, "geometry"] = route_linestring

        maritime_edges = gpd.GeoDataFrame(maritime_edges).set_crs(epsg=4326)

        maritime_nodes.to_parquet(output.nodes)
        maritime_edges.to_parquet(output.edges)


rule plot_maritime_network:
    input:
        nodes = "{OUTPUT_DIR}/maritime_network/nodes.gpq",
        edges = "{OUTPUT_DIR}/maritime_network/edges.gpq",
    output:
        edges_plot = "{OUTPUT_DIR}/maritime_network/edges.png",
    run:
        import geopandas as gpd
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np

        from open_gira.plot.utils import chop_at_antimeridian

        matplotlib.use("Agg")
        plt.style.use("bmh")

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world.geometry = world.geometry.boundary

        maritime_nodes = gpd.read_parquet(input.nodes)
        maritime_edges = gpd.read_parquet(input.edges)

        # whole network
        f, ax = plt.subplots(figsize=(16, 7))
        chop_at_antimeridian(maritime_edges, drop_null_geometry=True).plot(
            ax=ax,
            linewidth=0.5,
            alpha=0.8
        )
        world.plot(ax=ax, lw=0.5, alpha=0.2)
        ax.set_xticks(np.linspace(-180, 180, 13))
        ax.set_yticks([-60, -30, 0, 30, 60])
        ax.set_ylim(-65, 85)
        ax.set_xlim(-180, 180)
        ax.grid(alpha=0.3)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        f.savefig(output.edges_plot)


rule plot_port_connections:
    input:
        nodes = "{OUTPUT_DIR}/maritime_network/nodes.gpq",
        edges = "{OUTPUT_DIR}/maritime_network/edges.gpq",
    output:
        port_trade_plots = directory("{OUTPUT_DIR}/maritime_network/port_trade_plots"),
    run:
        import os

        import geopandas as gpd
        import matplotlib
        import matplotlib.pyplot as plt
        from tqdm import tqdm

        from open_gira.plot.utils import chop_at_antimeridian

        matplotlib.use("Agg")
        plt.style.use("bmh")

        maritime_nodes = gpd.read_parquet(input.nodes)
        maritime_edges = gpd.read_parquet(input.edges)

        # disambiguate the global view and plot the routes from each port, one port at a time
        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world.geometry = world.geometry.boundary

        ports = maritime_edges.from_port.unique()
        os.makedirs(output.port_trade_plots)
        for port_id in tqdm(ports):
            filepath = os.path.join(output.port_trade_plots, f"{port_id}.png")
            f, ax = plt.subplots(figsize=(10,10))
            port = maritime_nodes[maritime_nodes.id == f"{port_id}_land"]
            routes = maritime_edges[maritime_edges.from_port == port_id]
            if all(routes.geometry.isna()):
                continue
            try:
                chop_at_antimeridian(routes, drop_null_geometry=True).plot(
                    column="to_port",
                    categorical=True,
                    ax=ax
                )
            except ValueError:
                print(f"Failed to plot {port_id}, skipping...")
                continue
            maritime_nodes[maritime_nodes.id == f"{port_id}_land"].plot(
                ax=ax,
                markersize=500,
                marker="*",
                facecolor="none",
                color="r"
            )
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            world.plot(ax=ax, linewidth=0.5, alpha=0.4)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            port_name, = port.name
            ax.set_title(f"{port_id} ({port_name.replace('_', ', ')}) estimated routes")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            
            f.savefig(filepath)
            plt.close(f)