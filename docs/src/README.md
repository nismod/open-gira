# open-gira

This open-source [snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow can be used to analyse environmental risks to infrastructure
networks using global open data.

We can use `open-gira` to analyse open data on roads, railways and
their exposure to river flooding, or electricity transmission lines and how
they're affected by tropical cyclones.

`open-gira` is a work in progress.

## Goals
- Automated pipeline for reproducible analysis anywhere in the world
- Maps per-country and of larger areas
- Charts/stats of exposure per admin region, per hazard type, scenario, epoch
- Consider transport, electricity, water, communications systems
- Consider river flooding, storm surge coastal flooding, tropical cyclones
- Estimate direct damages to physical networks
- Estimate indirect effects of disruption - people affected, economic activity disrupted

## Non-goals
- Using closed data, which may be appropriate for other projects or use-cases
- Detailed operational/engineering level simulation
- Long-term planning
