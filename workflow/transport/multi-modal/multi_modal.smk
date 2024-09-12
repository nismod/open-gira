rule create_multi_modal_network:
    """
    Take previously created road, rail and maritime networks and combine them
    into a single multi-modal network with intermodal connections within
    distance limit of: any road node, any rail station and any maritime port.
    """
    input:
        admin_boundaries = "{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
        road_network_nodes = "{OUTPUT_DIR}/input/networks/road/{PROJECT_SLUG}/nodes.gpq",
        road_network_edges = "{OUTPUT_DIR}/input/networks/road/{PROJECT_SLUG}/edges.gpq",
        rail_network_nodes = "{OUTPUT_DIR}/input/networks/rail/{PROJECT_SLUG}/nodes.gpq",
        rail_network_edges = "{OUTPUT_DIR}/input/networks/rail/{PROJECT_SLUG}/edges.gpq",
        maritime_nodes = "{OUTPUT_DIR}/maritime_network/nodes.gpq",
        maritime_edges = "{OUTPUT_DIR}/maritime_network/edges.gpq",
    output:
        border_crossing_plot = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/border_crossings.png",
        nodes = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/nodes.gpq",
        edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/edges.gpq",
    script:
        "./multi_modal.py"


rule remove_edges_in_excess_of_threshold:
    """
    Take part of multi-modal network and remove edges that experience hazard
    values in excess of a given threshold.
    """
    input:
        edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/edges.gpq",
        raster = "{OUTPUT_DIR}/hazard/{HAZARD_SLUG}.tif",
    output:
        edges = temp("{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/{HAZARD_SLUG}/chunks/{CHUNK_SLUG}/edges.gpq")
    run:
        import geopandas as gpd
        import numpy as np

        from open_gira.disruption import filter_edges_by_raster

        i = int(wildcards.CHUNK_SLUG.split("-")[-1])
        edges = gpd.read_parquet(input.edges)
        edges = edges.loc[edges["mode"].isin({"road", "rail"}), :]
        chunk_size = int(np.ceil(len(edges) / workflow.cores))
        print(f"Chunk {i}: intersecting edges [{i * chunk_size}: {(i + 1) * chunk_size}]")

        surviving_edges = filter_edges_by_raster(
            edges.iloc[i * chunk_size: (i + 1) * chunk_size, :],
            input.raster,
            float(config["edge_failure_threshold"])
        )

        surviving_edges.to_parquet(output.edges)


rule join_intersection_results:
    """
    Pull together chunks of intersection operation to remove edges exceeding threshold value.
    """
    input:
        all_edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/edges.gpq",
        intersected_edges = expand(
            "{{OUTPUT_DIR}}/multi-modal_network/{{PROJECT_SLUG}}/{{HAZARD_SLUG}}/chunks/chunk-{chunk}/edges.gpq",
            chunk=range(0, workflow.cores)
        )
    output:
        edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/{HAZARD_SLUG}/edges.gpq",
    run:
        import geopandas as gpd
        import pandas as pd

        all_edges = gpd.read_parquet(input.all_edges)

        intersected_edges_post_hazard = pd.concat([gpd.read_parquet(path) for path in input.intersected_edges])
        edges = pd.concat([all_edges.loc[~all_edges["mode"].isin({"road", "rail"}), :], intersected_edges_post_hazard])

        edges.to_parquet(output.edges)