# Flow allocation

This workflow can allocate transport flows from origin-destination matrices over
networks. These networks may be intact or degraded.

## Description

1. We require an input trade origin-destination matrix (OD).

This should have the following format:
```
                   id partner_GID_0    value_kusd   volume_tons
thailand-latest_14_10           ABW  4.536817e-07  5.172909e-07
thailand-latest_14_10           AFG  3.524201e-06  8.233424e-07
thailand-latest_14_10           AGO  4.582354e-04  1.108041e-03
```

Where:
- `id` are node ids in the network to be routed over
- `partner_GID_0` is a partner country ISO 3-letter code
- `value_kusd` is the trade value in thousands of USD per year
- `volume_tons` is the trade mass (a.k.a. volume) in tonnes per year

This file should be located here:
```python
od = "{OUTPUT_DIR}/input/trade_matrix/{PROJECT_SLUG}/trade_nodes_total.parquet",
```

1. This workflow can route over intact or disrupted networks. If disrupting a
network, the raster used to represent the hazard must be available at:
```python
raster = "{OUTPUT_DIR}/hazard/{HAZARD_SLUG}.tif",
```

## Configuration

Review and amend `config/config.yaml`:
- `minimum_flow_volume_t`: discard flows in the OD with a `volume_tons` value less than this
- `edge_failure_threshold`: when failing a network subject to a hazard raster,
remove edges experiencing intensities greater than this


## Creation

Here's an example creation command where `PROJECT_SLUG` is `project-thailand`:
```bash
snakemake --cores all results/flow_allocation/project-thailand/routes_with_costs.pq",
```

And for the disrupted case, where `HAZARD_SLUG` is `hazard-thai-floods-2011`:
```bash
snakemake --cores all results/flow_allocation/project-thailand/hazard-thai-floods-2011/routes_with_costs.pq",
```
