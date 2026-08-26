"""
Network science type analysis of gridfinder electricity networks.
"""

rule calculate_network_metrics:
    """
    Find network attributes of gridfinder electricity networks. For example,
    edge and node country, betweenness centrality, average shortest path
    distance, assortativity, etc.

    Test with:
    snakemake -c1 results/power/by_country/PRI/network/attributes.pq
    """
    input:
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/nodes.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json",
    threads:
        16
    output:
        metrics="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/metrics.pq",
    run:
        import geopandas as gpd

        from open_gira.network_analysis import compute_all_metrics

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        edges = gpd.read_parquet(input.edges).set_crs(epsg=4326)
        nodes = gpd.read_parquet(input.nodes).set_crs(epsg=4326)

        metrics = compute_all_metrics(edges, nodes, n_workers=threads)
        metrics["iso_a3"] = wildcards.COUNTRY_ISO_A3
        metrics.to_frame().T.set_index("iso_a3").to_parquet(output.metrics)
