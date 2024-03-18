rule network_components:
    input:
        nodes="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/nodes.gpq",
        edges="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/edges.gpq",
    output:
        component_population="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/component_population.svg",
        component_map="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/network_map_by_component.png",
        component_data="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/components.parquet"
    script:
        "./network_components.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-road/component_population.svg
"""
