# Take split .geoparquet files for nodes or edges and output unified .geoparquet
# files for transport networks

rule join_network_nodes:
    input:
        lambda wildcards: expand(
            os.path.join("{OUTPUT_DIR}", "geoparquet", "{DATASET}_{FILTER_SLUG}", "processed", "slice-{i}_nodes_annotated.geoparquet"),
            **wildcards,
            i=range(config['slice_count'])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/nodes.geoparquet",
    script:
        "../scripts/join_data.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/nodes.geoparquet
"""


rule join_network_edges:
    input:
        # without a slice count in the output (we're aggregating), snakemake can't determine the slice in the input
        # therefore, use expand to generate the inputs from a list of slice numbers
        lambda wildcards: expand(
            os.path.join("{OUTPUT_DIR}", "geoparquet", "{DATASET}_{FILTER_SLUG}", "processed", "slice-{i}_edges_annotated.geoparquet"),
            **wildcards,
            i=range(config['slice_count'])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/edges.geoparquet"
    script:
        "../scripts/join_edges.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/edges.geoparquet
"""


rule assess_network_connectedness:
    input:
        nodes="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/nodes.geoparquet",
        edges="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/edges.geoparquet",
    output:
        component_population="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/component_population.pdf",
        component_map="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/network_map_by_component.png",
        component_data="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/components.parquet"
    script:
        "../scripts/assess_network_connectedness.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/component_population.pdf
"""
