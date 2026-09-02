"""
Network science type analysis of gridfinder electricity networks.
"""

from open_gira.io import cached_json_file_read


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


def network_metrics_affected_by_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their network metrics files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/metrics.pq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=country_set,  # list of str
    )


rule calculate_network_metrics_storm_set:
    """
    Concatenate network metrics for all countries potentially impacted by a
    storm set into a single table.

    Test with:
    snakemake -c1 results/power/by_storm_set/IBTrACS/network_metrics.pq
    """
    input:
        metrics = network_metrics_affected_by_storm_set
    output:
        network_attributes = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/network_metrics.pq",
    run:
        import pandas as pd

        metrics = pd.concat((pd.read_parquet(path) for path in input.metrics))
        metrics.to_parquet(output.network_attributes)
