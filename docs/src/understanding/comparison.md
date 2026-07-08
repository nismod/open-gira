# Which tool do I need?

open-gira is one of several open tools for climate and disaster risk to
infrastructure. They overlap, but they were built for different questions, and
the honest answer to "which should I use?" is often "more than one." This page
describes what open-gira is good for, what the neighbouring tools are good for,
and where they fit together.

> **A note on fairness.** The descriptions below reflect our understanding of
> each tool from its public documentation at the time of writing. Tools change,
> and we may have misread one. If you maintain a tool described here and
> something is wrong or out of date, please
> [open an issue](https://github.com/nismod/open-gira/issues) — we would rather
> correct it than mislead a reader.

## What open-gira is for

open-gira's distinctive feature is that it runs *end to end* from raw open data
to a risk estimate, as a single reproducible pipeline. You give it a place, a
network type, and a hazard; it downloads the source data, builds the network,
performs the hazard intersection, and produces damages and disruption — with the
dependency graph, the intermediate files, and the software environment all
pinned. There is no manual step where you prepare an exposure layer by hand and
feed it in.

This makes open-gira a good fit when you want to:

- go from **open source data to result** without assembling the exposure and
  network layers yourself;
- run the **same analysis over many places, hazards, scenarios, and epochs**,
  and have the bookkeeping handled by the workflow rather than by you;
- study **networked** infrastructure — where the question is not only "which
  assets are exposed?" but "what is disconnected when they fail?";
- reproduce or extend a published analysis exactly, from a lockfile and a
  configuration.

It is a less natural fit when you already have your own exposure and hazard
layers and simply want a damage calculation, or when you need a single
well-documented function to drop into an existing model. For those, the tools
below are often the better starting point.

## Neighbouring tools

### CLIMADA

[CLIMADA](https://climada.ethz.ch/) (ETH Zürich) is a mature, widely used
probabilistic risk-assessment platform built around the triad of *hazard*,
*exposure*, and *vulnerability*. It is strong on probabilistic event sets,
economic impact, and adaptation-option appraisal, and it has an extensive
scientific literature behind it.

CLIMADA and open-gira meet at the exposure layer. CLIMADA typically represents
exposure as points or gridded values; open-gira builds explicit network
topology and can reason about connectivity and cascading service loss. A natural
division of labour is to use open-gira to construct networks and hazard
intersections and CLIMADA for probabilistic impact and adaptation appraisal. We
regard interoperability with CLIMADA (for example, exporting open-gira exposure
as CLIMADA `Exposures`) as a desirable direction rather than a competitive one.

### RA2CE

[RA2CE](https://github.com/Deltares/ra2ce) (Deltares) is a resilience assessment
tool focused on transport networks and criticality — origin–destination routing,
redundancy, and the loss of accessibility when links fail. Its focus on
transport-network criticality is close to open-gira's interest in disruption,
and the two make a sensible pairing: RA2CE has a well-developed accessibility
and detour analysis; open-gira has an automated multi-hazard, multi-region data
pipeline. RA2CE has a track record of country-scale deployments, which is a
strength open-gira does not yet match.

### physrisk

[physrisk](https://github.com/os-climate/physrisk) (OS-Climate) is a calculation
engine for the physical climate risk of financial assets and portfolios. Its
audience is largely financial — asset-level damage and business interruption
aggregated to a portfolio. open-gira's audience is largely infrastructure and
development analysis, and its unit of interest is the network rather than the
balance sheet. The two are more complementary than overlapping; shared hazard
data and metadata conventions (see below) are the most likely point of contact.

### DamageScanner

[DamageScanner](https://github.com/VU-IVM/DamageScanner) is a focused,
lightweight Python library for object- and raster-based damage assessment given
hazard, exposure, and vulnerability curves. It does one part of the problem —
the damage calculation — cleanly and with little ceremony. If you already have
your layers and want a transparent damage function, DamageScanner is a smaller
and simpler dependency than open-gira. open-gira solves a larger problem (data
acquisition, network construction, orchestration) at correspondingly greater
weight.

## A rough guide

| If you want to…                                                             | Consider           |
| --------------------------------------------------------------------------- | ------------------ |
| Go from open data to network risk, reproducibly, over many places           | **open-gira**      |
| Do probabilistic impact and adaptation appraisal on exposure you provide     | CLIMADA            |
| Assess transport-network criticality and accessibility loss                  | RA2CE              |
| Compute physical climate risk for a financial portfolio                      | physrisk           |
| Apply a damage curve to layers you already have                              | DamageScanner      |

These are starting points, not walls. open-gira is designed to sit within this
ecosystem — consuming shared hazard data, and (as a development direction)
exporting to the formats these tools expect — rather than to replace any one of
them.

## Interoperability

open-gira produces standard geospatial outputs (GeoParquet, GeoPackage, raster)
in documented coordinate reference systems, which most of the tools above can
read. Two conventions matter for interchange:

- **[Risk Data Library Standard (RDLS)](https://docs.riskdatalibrary.org/)**
  metadata for describing hazard and risk data, so that open-gira outputs are
  self-describing to other tools and catalogues.
- **CLIMADA `Exposures`** as an export target, so open-gira networks can be
  carried into CLIMADA's probabilistic machinery.

Where these are implemented, they are documented in
[Using open-gira data](data/sources-and-licences.md); where they are planned
rather than present, we say so there.
