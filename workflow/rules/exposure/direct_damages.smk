# estimate damage fraction for exposed assets
# damage fraction is a function of hazard intensity (expressed as damage curves)


rule direct_damages:
    input:
        exposure = "{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/slice-{i}.parquet",
    output:
        damages = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/slice-{i}.parquet",
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        # TODO: determine the hazard type from the hazard slug?
        hazard_type="flood",
    script:
        "../../scripts/transport/direct_damages.py"


"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-5.parquet
"""


rule plot_damage_fractions:
    input:
        damages = rules.join_direct_damages.output.joined
    output:
        plots = directory("{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/plots")
    script:
        "../../scripts/transport/plot_damage_fractions.py"


"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/plots
"""
