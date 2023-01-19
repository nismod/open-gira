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
        gridfinder="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
        targets="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        powerplants="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
    run:
        import os

        import geopandas as gpd

        gridfinder = gpd.read_parquet(input.gridfinder)
        targets = gpd.read_parquet(input.targets)
        powerplants = gpd.read_parquet(input.powerplants)

        countries = gpd.read_parquet(input.admin_bounds).rename(columns={"GID_0": "iso_a3"})
        country = countries[countries.iso_a3 == wildcards.COUNTRY_ISO_A3]
        country = country[["iso_a3", "geometry"]]

        os.makedirs(os.path.dirname(output.gridfinder), exist_ok=True)

        country_grid = gridfinder.sjoin(country, how="inner")
        country_grid[gridfinder.columns].to_parquet(output.gridfinder)
        country_targets = targets.sjoin(country, how="inner")
        country_targets[targets.columns].to_parquet(output.targets)
        country_powerplants = powerplants.sjoin(country, how="inner")
        country_powerplants[powerplants.columns].to_parquet(output.powerplants)

"""
Test with:
snakemake -c1 results/power/country/HTI/grid.geoparquet
"""


rule create_power_network:
    """
    Combine power plant, consumer and transmission data for given area
    """
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
        targets="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        gridfinder="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
    output:
        edges="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/country/{COUNTRY_ISO_A3}/network/nodes.geoparquet",
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
