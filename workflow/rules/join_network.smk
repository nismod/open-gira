# Take split .geoparquet files for nodes or edges and output unified .geoparquet
# files for transport networks

rule join_network_nodes:
    input:
        lambda wildcards: expand(
            os.path.join("{OUTPUT_DIR}", "geoparquet", "{DATASET}_{FILTER_SLUG}", "slice-{i}_road_nodes_annotated.geoparquet"),
            **wildcards,
            i=range(config['slice_count'])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_nodes.geoparquet",
    script:
        "../scripts/join_data.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/road_nodes.geoparquet
"""


rule join_network_edges:
    input:
        # without a slice count in the output (we're aggregating), snakemake can't determine the slice in the input
        # therefore, use expand to generate the inputs from a list of slice numbers
        lambda wildcards: expand(
            os.path.join("{OUTPUT_DIR}", "geoparquet", "{DATASET}_{FILTER_SLUG}", "slice-{i}_road_edges_annotated.geoparquet"),
            **wildcards,
            i=range(config['slice_count'])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_edges.geoparquet"
    script:
        "../scripts/join_edges.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/road_edges.geoparquet
"""


rule plot_network_connectedness:
    input:
        nodes="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_nodes.geoparquet",
        edges="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_edges.geoparquet",
    output:
        component_population="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_component_population.pdf",
        component_map="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/road_network_map_by_component.png"
    script:
        "../scripts/plot_network_connectedness.py"

"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core/road_component_population.pdf
"""
