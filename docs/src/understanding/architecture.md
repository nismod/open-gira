# Architecture overview

This page describes how open-gira is put together: the shape of the workflow,
the sectors and hazards it covers, the outputs it leaves on disk, and the
libraries it rests on. It is meant to give a reader enough of a mental model to
navigate the rest of the documentation and the code.

## The workflow is a graph of rules

open-gira is a [snakemake](https://snakemake.readthedocs.io/) workflow. Rather
than a program you run start-to-finish, it is a collection of *rules*, each of
which knows how to turn some input files into some output files. You ask for a
result file; snakemake works backwards to determine which rules must run, in
what order, to produce it, and then runs them.

The consequences of this design are worth stating plainly, because they shape
everything else:

- **The interface is the set of output filenames.** You do not call functions;
  you request paths. `results/wales-latest_filter-road-primary/edges.gpq` is a
  request for the edges of the primary road network in Wales. The meaning is
  encoded in the path. This is powerful and terse, and it is also the main thing
  a new user has to learn — see [Targets and wildcards](../reference/targets-and-wildcards.md).
- **Work is only done once.** If an intermediate file already exists and is up
  to date, its rule does not re-run. Analyses are incremental by construction.
- **The dependency graph is explicit.** You can ask snakemake to show what it
  *would* do (`snakemake -n`) or to draw the graph, before committing to a long
  computation.

The rules are defined across roughly sixty `.smk` files under `workflow/`, wired
together by `workflow/Snakefile`. They are grouped by sector and by the coupling
of a sector to a hazard.

## The sector × hazard matrix

open-gira's rule files fall into three groups: those that prepare a **sector**
(building networks and exposure), those that prepare a **hazard**, and those
that **couple** a sector to a hazard to produce risk. Reading `workflow/Snakefile`
top to bottom is reading exactly this structure.

**Sectors** (network and exposure construction):

- `transport/` — roads and railways from OpenStreetMap, plus maritime and
  multi-modal network assembly and trade-flow routing;
- `power/` — electricity transmission networks from
  [gridfinder](https://github.com/carderne/gridfinder) and power plants from WRI;
- `buildings/` — built-area and building-exposure value from GHSL, GIRI, and
  economic capital-stock sources;
- `nature-ecosystems/` — land cover, hydrobasins, and nature-based-solution
  suitability;
- `population-economy/` — gridded population and GDP.

**Hazards**:

- `flood/` — river flooding (JRC, Aqueduct) and coastal flooding (Deltares);
- `tropical-cyclone/` — wind fields derived from several track sets (IBTrACS,
  STORM, IRIS, CHAZ, and Emanuel tracks);
- `landslide/` — landslide susceptibility (Arup).

**Couplings** (sector × hazard → risk):

| Sector ↓ / Hazard →   | Flood                 | Tropical cyclone    | Landslide               |
| --------------------- | --------------------- | ------------------- | ----------------------- |
| Transport             | `transport-flood/`    | —                   | `transport-landslide/`  |
| Power (electricity)   | `power-flood/`        | `power-tc/`         | —                       |

The blank cells are honest: not every sector is coupled to every hazard. The
most developed couplings are transport/flooding and electricity/tropical-cyclone,
which are also the two [tutorials](../tutorials/wales-roads-flooding.md) and the
two [risk-analysis usage pages](../user-guide/usage/risk-analysis.md).

A coupling rule generally does three things: intersect the network with the
hazard (which network elements experience what intensity), apply a
vulnerability or damage curve (intensity → damage fraction → cost), and
aggregate (to administrative regions, return periods, or expected annual
damage).

## What ends up on disk

All outputs live under a directory whose name ends in `results` (this suffix is
enforced; see [constraints and gotchas](../reference/constraints-and-gotchas.md)).
Within it, paths encode the analysis. A representative slice:

```
results/
├── <dataset>_<filter-slug>/          # a network, e.g. wales-latest_filter-road-primary
│   ├── edges.gpq                      #   network edges (GeoParquet)
│   └── nodes.gpq                      #   network nodes
├── input/                             # downloaded and pre-processed source data
├── power/                             # electricity network and derived products
├── exposure/                          # network × hazard intersections
├── direct_damages/                    # damage fractions and costs
└── ...                                # further products, per workflow
```

The exact tree depends on which targets you request; the workflow only creates
what is needed. The [targets and wildcards reference](../reference/targets-and-wildcards.md)
is the map from a rule to the paths it produces.

## The library stack

open-gira is deliberately thin over a few focused libraries, several of which
were developed alongside it:

- **[snail](https://github.com/nismod/snail)** — vector–raster intersection: it
  answers "which network segments fall in which hazard raster cells, at what
  intensity?" Used throughout the hazard-intersection rules.
- **[snkit](https://github.com/tomalrussell/snkit)** — network cleaning and
  assembly: turning raw geometries into topologically sound node/edge networks.
- **[irv-datapkg](https://pypi.org/project/irv-datapkg/)** — packaging and
  cropping of gridded source datasets to administrative boundaries.

Around these sit the general geospatial and scientific Python stack (GeoPandas,
rasterio, xarray, GDAL, and others), pinned via [pixi](https://pixi.sh/) for
reproducibility. The full list is in `pixi.toml`.

Keeping domain logic in `snail` and `snkit` — libraries with their own tests and
releases — means open-gira itself is mostly orchestration: configuration
checking, rule definitions, and the scripts that glue library calls together.
