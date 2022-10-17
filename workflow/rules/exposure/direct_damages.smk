# estimate damage fraction for exposed assets
# damage fraction is a function of hazard intensity (expressed as damage curves)


rule direct_damages:
    input:
        exposure = rules.network_raster.output.geoparquet
    output:
        damage_fraction = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/fraction/{SLICE_SLUG}.geoparquet",
        damage_cost = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/cost/{SLICE_SLUG}.geoparquet",
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        # determine the hazard type from the hazard slug, e.g. flood, earthquake, storm
        hazard_type=lambda wildcards: config["hazard_types"][wildcards.HAZARD_SLUG.replace('hazard-', '')]
    script:
        "../../scripts/transport/direct_damages.py"


"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/cost/slice-5.geoparquet
"""


rule plot_damage_distributions:
    input:
        damages = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction.geoparquet"
    output:
        plots = directory("{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/fraction_plots")
    script:
        "../../scripts/transport/plot_damage_distributions.py"


"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/fraction_plots
"""
