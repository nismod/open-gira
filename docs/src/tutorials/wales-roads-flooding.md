# Flood risk to Wales' roads

This tutorial takes you from a clean checkout to a first result: expected annual
flood damages to the primary road network of Wales, aggregated to the country.
It is meant to run on a laptop, and every command is meant to be copied and
pasted. If a step does something surprising, that is a documentation bug — please
[tell us](https://github.com/nismod/open-gira/issues).

Wales is a good first target because it is small: the road network builds in a
short time, and the intermediate files are modest. Once this works, the same
commands work for anywhere else by changing the place in the target path.

> **On timings and sizes.** Where this page gives a runtime or a disk figure, it
> is a rough guide on a typical laptop, not a guarantee — hazard downloads in
> particular dominate, and depend on your connection. Use the dry-run trick in
> each step to see what *would* happen before committing to it.
>
> *(maintainer to-do: replace the qualitative figures below with measured
> values from a reference machine, and pin the exact hazard dataset used.)*

## Before you start

You will need open-gira installed. If you have not done that yet, follow
[Installation](../user-guide/installation.md); the short version, using
[pixi](https://pixi.sh/), is:

```bash
git clone https://github.com/nismod/open-gira.git
cd open-gira
pixi install
pixi shell        # gives you a shell with the environment active
```

Every command below assumes you are inside `pixi shell`, in the repository root.

## Step 1 — Look before you leap

open-gira never makes you guess what a command will do. Before running anything,
ask snakemake to *plan* the work with `-n` (a dry run):

```bash
snakemake -n --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

This prints the rules that would run to produce the Welsh primary-road network,
without running them. Read the path like a sentence: `wales-latest` is the
OpenStreetMap extract, and `filter-road-primary` selects the primary-road filter
defined in `config/osm_filters/road-primary.txt`. The
[targets and wildcards reference](../reference/targets-and-wildcards.md) explains
this grammar in full.

## Step 2 — Build the network

Now run it for real:

```bash
snakemake --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

snakemake downloads the Wales OpenStreetMap extract, filters it to primary
roads, slices it, converts it to GeoParquet, cleans and annotates it, and joins
the pieces back together (the [road page](../user-guide/usage/network-creation/road.md)
describes each step). This should take a few minutes on a laptop.

When it finishes you have two files:

```
results/wales-latest_filter-road-primary/
├── edges.gpq     # the road segments
└── nodes.gpq     # the junctions
```

You can open these directly in QGIS, or in Python:

```python
import geopandas as gpd

edges = gpd.read_parquet("results/wales-latest_filter-road-primary/edges.gpq")
print(len(edges), "edges")
edges.plot()
```

## Step 3 — Choose a flood hazard

Damage estimation needs a hazard. open-gira reads hazard datasets from
`config/config.yaml` under `hazard_datasets`, with a matching entry in
`hazard_types`. The [transport/flooding page](../user-guide/usage/risk-analysis/transport-flooding.md)
describes the configuration in detail.

For this tutorial, use a river-flood dataset (for example an Aqueduct river set).
Confirm it is configured:

```bash
# there should be a hazard-<name> entry you can request; check the config
grep -A3 hazard_datasets config/config.yaml
```

> **Be aware of download size.** Global flood-hazard raster sets can be large —
> potentially several gigabytes before cropping. open-gira crops to the extent
> of your network, but it downloads the source first. Run the dry run in the
> next step to see what will be fetched.

## Step 4 — Compute expected annual damages

Ask for the damages aggregated to country level (administrative level 0). Do the
dry run first:

```bash
snakemake -n --cores 2 -- \
  results/wales-latest_filter-road-primary/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
```

Read the target path as the chain of the analysis: the Welsh primary-road
network, intersected with the `aqueduct-river` hazard, reduced to expected
annual damages and cost per return period, summed, and aggregated to admin level
0. When the plan looks right, drop the `-n`:

```bash
snakemake --cores 2 -- \
  results/wales-latest_filter-road-primary/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
```

Behind that single request, open-gira crops the hazard to Wales, splits the road
edges on the raster grid so no edge spans two pixels, reads the flood depth for
each split edge at each return period, applies a damage curve per asset type,
converts damage fractions to cost using rehabilitation-cost estimates, and
integrates over return periods to get expected annual damages.

## Step 5 — Look at the result

```python
import geopandas as gpd

ead = gpd.read_parquet(
    "results/wales-latest_filter-road-primary/hazard-aqueduct-river/"
    "EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet"
)
print(ead.columns.tolist())
ead.plot(column="EAD", legend=True)
```

You now have an expected-annual-damage figure for Welsh primary roads under the
chosen flood hazard, and a map to go with it.

## Where to go next

- Change `wales-latest` to another entry from `infrastructure_datasets` in
  `config/config.yaml` to run the same analysis elsewhere.
- Change `admin-level-0` to `admin-level-1` for a finer regional breakdown.
- Read [Targets and wildcards](../reference/targets-and-wildcards.md) to see the
  full space of things you can request.
- Try the second tutorial,
  [Cyclone risk to Puerto Rico's grid](puerto-rico-grid-cyclone.md), for the
  electricity/tropical-cyclone workflow.
