
rule landslide_direct_damages:
    input:
        unsplit = rules.create_transport_network.output.edges,  # for pre-intersection geometry
        exposure = rules.rasterise_osm_network.output.geoparquet,
        rehab_cost=lambda wildcards: f"config/rehab_costs/{wildcards.FILTER_SLUG.split('-')[1]}.csv",
        damage_curves="config/damage_curves/",
    output:
        split_ead_and_cost_per_trigger = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/split_EAD_and_cost_per_trigger/{SLICE_SLUG}.geoparquet",
        ead_and_cost_per_trigger = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_trigger/{SLICE_SLUG}.geoparquet",
    params:
        # determine the network type from the filter, e.g. filter-road-primary: road, filter-rail: rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.split("-")[1],
    script:
        "./landslide_direct_damages.py"
