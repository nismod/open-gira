rule allocate_intact_network:
    """
    Allocate a trade OD matrix across a multi-modal transport network

    Test with:
    snakemake -c1 -- results/flow_allocation/project-thailand/edges.gpq 
    """
    input:
        edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/edges.gpq",
        od = "{OUTPUT_DIR}/input/trade_matrix/{PROJECT_SLUG}/trade_nodes_total.parquet",
    threads: workflow.cores
    params:
        # if this changes, we want to trigger a re-run
        minimum_flow_volume_t = config["minimum_flow_volume_t"],
    output:
        routes = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/routes.pq",
        edges_with_flows = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/edges.gpq",
    script:
        "./allocate.py"


rule allocate_degraded_network:
    """
    Allocate a trade OD matrix across a multi-modal transport network which has
    lost edges as a result of intersection with a hazard map.

    Test with:
    snakemake -c1 -- results/flow_allocation/project-thailand/hazard-thai-floods-2011-JBA/edges.gpq 
    """
    input:
        edges = "{OUTPUT_DIR}/multi-modal_network/{PROJECT_SLUG}/{HAZARD_SLUG}/edges.gpq",
        od = "{OUTPUT_DIR}/input/trade_matrix/{PROJECT_SLUG}/trade_nodes_total.parquet",
    threads: workflow.cores
    params:
        # if this changes, we want to trigger a re-run
        minimum_flow_volume_t = config["minimum_flow_volume_t"],
    output:
        routes = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/{HAZARD_SLUG}/routes.pq",
        edges_with_flows = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/{HAZARD_SLUG}/edges.gpq",
    script:
        "./allocate.py"


rule accumulate_route_costs_intact:
    """
    For each route in the OD (source -> destination pair), lookup the edges of
    the least cost route (the route taken) and sum those costs. Store alongside
    value and volume of route.
    
    Test with:
    snakemake -c1 -- results/flow_allocation/project-thailand/routes_with_costs.pq
    """
    input:
        routes = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/routes.pq",
        edges_with_flows = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/edges.gpq",
    output:
        routes_with_costs = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/routes_with_costs.pq",
    run:
        from open_gira.routing import lookup_route_costs

        lookup_route_costs(input.routes, input.edges_with_flows).to_parquet(output.routes_with_costs)


rule accumulate_route_costs_degraded:
    """
    For each route in the OD (source -> destination pair), lookup the edges of
    the least cost route (the route taken) and sum those costs. Store alongside
    value and volume of route.
    
    Test with:
    snakemake -c1 -- results/flow_allocation/project-thailand/hazard-thai-floods-2011-JBA/routes_with_costs.pq
    """
    input:
        routes = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/{HAZARD_SLUG}/routes.pq",
        edges_with_flows = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/{HAZARD_SLUG}/edges.gpq",
    output:
        routes_with_costs = "{OUTPUT_DIR}/flow_allocation/{PROJECT_SLUG}/{HAZARD_SLUG}/routes_with_costs.pq",
    run:
        from open_gira.routing import lookup_route_costs 

        lookup_route_costs(input.routes, input.edges_with_flows).to_parquet(output.routes_with_costs)
