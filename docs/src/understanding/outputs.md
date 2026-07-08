# What open-gira produces

This page is for the reader deciding whether open-gira is worth their time: what
it produces, who has used it, and what evidence there is that it works. We try to
show rather than assert, and to distinguish what we can point to from what is
still in progress.

> **Draft note for maintainers.** Several items below are marked *(to add)*.
> These are real outputs and deployments that belong on this page; they need a
> specific link, citation, or figure before publication. Please fill them in
> rather than leaving the placeholder.

## Outputs

A single open-gira run produces some combination of:

- **Networks** — cleaned, topologically connected road, rail, or electricity
  networks as GeoParquet (`nodes.gpq`, `edges.gpq`), ready to load in GeoPandas
  or QGIS.
- **Exposure** — tables recording which network elements intersect which
  hazard, at what intensity, per return period and scenario.
- **Direct damages** — damage fractions and reconstruction costs per asset,
  and expected annual damages (EAD) aggregated to administrative regions.
- **Disruption** — for the electricity/cyclone workflow, estimates of customers
  or population disconnected, by event or return period.
- **Maps and charts** — rendered figures of the above, for inspection and
  reporting.

The [tutorials](../tutorials/wales-roads-flooding.md) produce concrete examples
of each on a laptop.

## Releases and reproducibility

open-gira has been released regularly since 2023. The releases are not
cosmetic — each corresponds to real capability, and the record is public on the
[releases page](https://github.com/nismod/open-gira/releases):

| Version        | Date       | Notable content                                             |
| -------------- | ---------- | ----------------------------------------------------------- |
| v0.1.0         | 2023-09    | Initial release: road/flood, rail/flood, cyclone/grid       |
| v0.2.0         | 2024-04    | Event-set flood damages; composite networks; grid disruption |
| v0.3.0         | 2024-09    | Multi-modal (road/rail/maritime) network integration        |
| v0.3.1         | 2024-12    | Trade-flow routing; buildings exposure; STORM return periods |
| v0.3.2         | 2025-03    | Nature-based-solution opportunity areas (GCA-funded)        |
| v0.4.0-alpha   | 2025-03    | Power network / flood intersection; Apple-silicon support   |
| v0.4.1         | 2025-12    | Move to pixi lockfile; CHAZ cyclone tracks; faster disruption |
| v0.4.2         | 2026-03    | SLURM cluster execution                                      |

Two things about this record are worth drawing out for an evaluator:

- **The outputs are citable.** open-gira and its data releases carry a Zenodo
  DOI ([10.5281/zenodo.14537079](https://doi.org/10.5281/zenodo.14537079)). See
  [How to cite open-gira](../citation.md).
- **The pipeline is tested and pinned.** Every change runs an integration test
  suite that exercises the workflow rules against a small sample dataset, and
  the software environment is pinned with a [pixi](https://pixi.sh/) lockfile,
  so a given release can be reproduced rather than merely re-run.

## Where it has been used

open-gira has supported a number of studies and deployments. The following are
known to the maintainers; specifics are being collected here.

- **Global Center on Adaptation — nature-based solutions.** The nature-based
  solutions opportunity-area work (v0.3.2) was carried out as part of the GCA
  project *Scaling Investments in Nature-based Solutions for Climate-Resilient
  Infrastructure*.
- **Global infrastructure risk visualisation.** Outputs have fed the
  [Global Infrastructure Risk Model and Resilience Index](https://global.infrastructureresilience.org/)
  visualisation platform. *(confirm exact relationship and link)*
- **Published studies.** *(to add: citations of papers that used open-gira,
  with DOIs)*
- **Country and regional studies.** *(to add: named studies and the sectors /
  hazards they covered)*
- **External use.** open-gira has been forked and adapted for analyses beyond
  the core team. *(to add: examples, e.g. sovereign-risk analyses)*

## Usage metrics

To describe adoption honestly we prefer figures that reflect use over figures
that reflect popularity. The following are worth tracking, and we intend to
report them here as a baseline accrues:

- Zenodo download counts for the software and data releases;
- documentation-site analytics;
- citations of the *outputs* (data releases), not only the repository;
- an inventory of public forks and derived analyses.

*(to add: current figures, with the date they were measured)*
