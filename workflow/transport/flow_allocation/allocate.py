
import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from open_gira.routing import route_from_all_nodes, RouteResult


if __name__ == "__main__":

    print("Reading network...")
    # read in global multi-modal transport network
    edges = gpd.read_parquet(snakemake.input.edges)
    available_destinations = edges[edges["mode"] == "imaginary"].to_id.unique()
    available_country_destinations = [d.split("_")[-1] for d in available_destinations if d.startswith("GID_")]

    print("Reading OD matrix...")
    # read in trade OD matrix
    od = pd.read_parquet(snakemake.input.od)
    print(f"OD has {len(od):,d} flows")

    # 5t threshold drops THL road -> GID_0 OD from ~21M -> ~2M
    minimum_flow_volume_tons = snakemake.config["minimum_flow_volume_t"]
    od = od[od.volume_tons > minimum_flow_volume_tons]
    print(f"After dropping flows with volume < {minimum_flow_volume_tons}t, OD has {len(od):,d} flows")

    # drop any flows we can't find a route to
    od = od[od.partner_GID_0.isin(available_country_destinations)]
    print(f"After dropping unrouteable destination countries, OD has {len(od):,d} flows")

    routes: RouteResult = route_from_all_nodes(od, edges, snakemake.threads)

    print("Writing routes to disk as parquet...")
    pd.DataFrame(routes).T.to_parquet(snakemake.output.routes)

    print("Assigning route flows to edges...")
    edges["value_kusd"] = 0
    edges["volume_tons"] = 0
    value_col_id = edges.columns.get_loc("value_kusd")
    volume_col_id = edges.columns.get_loc("volume_tons")
    for (from_node, destination_country), route_data in tqdm(routes.items()):
        edges.iloc[route_data["edge_indices"], value_col_id] += route_data["value_kusd"]
        edges.iloc[route_data["edge_indices"], volume_col_id] += route_data["volume_tons"]

    print("Writing edge flows to disk as geoparquet...")
    edges.to_parquet(snakemake.output.edges_with_flows)

    print("Done")