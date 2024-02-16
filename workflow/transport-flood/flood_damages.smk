"""
Estimate damage fraction for exposed assets
Damage fraction is a function of hazard intensity (expressed as damage curves)
"""


rule return_period_direct_damages:
    input:
        unsplit = rules.create_transport_network.output.edges,  # for pre-intersection geometry
        exposure = rules.rasterise_osm_network.output.geoparquet,
        rehab_cost=lambda wildcards: f"config/rehab_costs/{wildcards.FILTER_SLUG.replace('filter-', '')}.csv",
        damage_curves="config/damage_curves/",
    output:
        damage_fraction = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/fraction_per_RP/{SLICE_SLUG}.geoparquet",
        damage_cost = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/cost_per_RP/{SLICE_SLUG}.geoparquet",
        expected_annual_damages = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD/{SLICE_SLUG}.geoparquet",
        return_period_and_ead = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_RP/{SLICE_SLUG}.geoparquet",
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        # determine the hazard type from the hazard slug, e.g. flood, earthquake, storm
        hazard_type=lambda wildcards: config["hazard_types"][wildcards.HAZARD_SLUG.replace('hazard-', '')]
    script:
        "./return_period_direct_damages.py"

"""
Test with:
snakemake --cores 1 results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/slice-5.geoparquet
"""


rule plot_damage_distributions:
    input:
        damages = "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction_per_RP.geoparquet"
    output:
        plots = directory("{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction_plots")
    script:
        "./plot_damage_distributions.py"

"""
Test with:
snakemake --cores 1 results/egypt-latest_filter-road/hazard-aqueduct-river/damage_fraction_plots
"""


rule event_set_direct_damages:
    input:
        unsplit = rules.create_transport_network.output.edges,  # for pre-intersection geometry
        exposure = rules.rasterise_osm_network.output.geoparquet,
        rehab_cost=lambda wildcards: f"config/rehab_costs/{wildcards.FILTER_SLUG.replace('filter-', '')}.csv",
        damage_curves="config/damage_curves/",
    output:
        damage_fraction = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/fraction/{SLICE_SLUG}.gpq",
        damage_cost = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/cost/{SLICE_SLUG}.gpq",
    params:
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', ''),
        hazard_type=lambda wildcards: config["hazard_types"][wildcards.HAZARD_SLUG.replace('hazard-', '')]
    script:
        "./event_set_direct_damages.py"

"""
Test with:
snakemake --cores 1 results/direct_damages/thailand-latest_filter-road/hazard-jba-event/cost/slice-5.gpq
"""


rule concat_event_set_direct_damages:
    input:
        slices = lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "direct_damages",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "{COST_OR_FRACTION}",
                "slice-{i}.gpq",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        damage_cost = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{COST_OR_FRACTION}.gpq",
    run:
        import geopandas as gpd
        import pandas as pd

        slices = []
        for path in input.slices:
            slices.append(gpd.read_parquet(path))
        pd.concat(slices).to_parquet(output.damage_cost)

"""
Test with:
snakemake -c1 results/direct_damages/thailand-latest_filter-road/hazard-jba-event/cost.gpq
"""