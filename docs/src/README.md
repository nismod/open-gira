# open-gira

open-gira is an open-source [snakemake](https://snakemake.readthedocs.io/)
workflow for analysing environmental risks to infrastructure networks using
global open data. Given a place and a hazard, it builds an infrastructure
network from open data, works out where the hazard intersects the network, and
estimates the resulting damage and disruption.

Concretely, open-gira can:

- build road, rail, and electricity networks from
  [OpenStreetMap](https://www.openstreetmap.org/) and other open sources;
- intersect those networks with hazard data for river flooding, coastal
  flooding, tropical cyclones, and landslides;
- estimate direct damages to the physical assets, and the disruption that
  follows (people disconnected, service lost);
- do this reproducibly, anywhere in the world, from raw data to result.

The name expands to *Open Global Infrastructure Risk/Resilience Analysis*.

## What this documentation is for

Readers arrive at open-gira for different reasons, and this book is organised
around them.

- If you are **deciding whether open-gira is the right tool**, start with
  [Which tool do I need?](understanding/comparison.md), which sets open-gira
  beside neighbouring tools, and [What open-gira produces](understanding/outputs.md),
  which shows the outputs and the work that has used them.
- If you want to **run an analysis**, the [tutorials](tutorials/wales-roads-flooding.md)
  walk from a clean checkout to a first result on a laptop. The
  [usage guide](user-guide/usage.md) then documents each workflow in turn.
- If you want to **use open-gira's outputs** without running the pipeline
  yourself, see [Data sources and licences](data/sources-and-licences.md).
- If you want to **extend or contribute**, see
  [Contributing to open-gira](contributing/index.md).

## Status

open-gira is research software under active development. It is versioned and
released, its outputs carry a citable DOI, and its workflows are exercised by
an automated test suite on every change (see the badges in the repository
[README](https://github.com/nismod/open-gira)). It has been used in a number of
published studies and deployments, described in
[What open-gira produces](understanding/outputs.md).

That said, coverage is uneven: some sectors and hazards are more mature than
others, and some documented features are ahead of the code or behind it. Where
we know of a limitation, we try to say so plainly, both in this book and in the
[constraints and gotchas](reference/constraints-and-gotchas.md) reference.

## Scope

open-gira aims to be:

- an automated, reproducible pipeline that runs anywhere in the world;
- a producer of maps, charts, and tables of exposure and risk, broken down by
  administrative region, hazard, scenario, and epoch;
- a tool spanning multiple systems (transport, electricity, and, increasingly,
  buildings and ecosystems) and multiple hazards.

It is deliberately *not*:

- a user of closed or proprietary data (appropriate for other projects, but not
  this one);
- an operational or engineering-level simulator;
- a long-term infrastructure planning tool.

## Acknowledgements

This research has received funding from the FCDO Climate Compatible Growth
Programme, the World Bank Group, the UK Natural Environment Research Council
(NERC) through the UK Centre for Greening Finance and Investment (CGFI), and the
Global Center on Adaptation (GCA). The views expressed here do not necessarily
reflect the funders' official policies.
