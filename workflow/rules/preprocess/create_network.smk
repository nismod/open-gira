"""
Network creation routines
"""


rule subset_grid_inputs_by_country:
    """
    Subset the gridfinder, targets and powerplants datasets to country bounds
    """
    input:
        gridfinder="{OUTPUT_DIR}/power/gridfinder.geoparquet",
        targets="{OUTPUT_DIR}/power/targets.geoparquet",
        powerplants="{OUTPUT_DIR}/power/powerplants.geoparquet",
        admin_bounds="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet"
    output:
        gridfinder="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
        targets="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        powerplants="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
    resources:
        mem_mb=8192
    script:
        "../../scripts/preprocess/slice_network_assets.py"

"""
Test with:
snakemake -c1 results/power/by_country/HTI/grid.geoparquet
"""


rule create_power_network:
    """
    Combine power plant, consumer and transmission data for given area
    """
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
        targets="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        gridfinder="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
    output:
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/nodes.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json",
    script:
        "../../scripts/preprocess/create_electricity_network.py"

"""
Test with:
snakemake -c1 results/power/edges.geoparquet
"""


rule create_transport_network:
    """
    Take .geoparquet OSM files and output files of cleaned network nodes and edges
    """
    conda: "../../../environment.yml"
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
        "../../scripts/transport/create_{params.network_type}_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-road/processed/slice-0_edges.geoparquet
"""
