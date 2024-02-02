

rule create_transport_network:
    """
    Take .geoparquet OSM files and output files of cleaned network nodes and edges
    """
    input:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_edges.geoparquet",
        admin="{OUTPUT_DIR}/input/admin-boundaries/gadm36_levels.gpkg",
    output:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_edges.geoparquet"
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        # pass in the slice number so we can label edges and nodes with their slice
        # edge and node IDs should be unique across all slices
        slice_number=lambda wildcards: int(wildcards.SLICE_SLUG.replace('slice-', ''))
    script:
        # template the path string with a value from params (can't execute .replace in `script` context)
        "./create_{params.network_type}_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-road/processed/slice-0_edges.geoparquet
"""
