# Targets and wildcards

In open-gira you do not call functions — you request files. The path you ask for
*is* the specification of the analysis. This page is the reference for that
grammar: the wildcards a path is built from, what each one may contain, and how
they compose into a target.

> **Keeping this page honest.** The wildcard definitions below are transcribed
> from the `wildcard_constraints` block in `workflow/Snakefile`. Hand-copied,
> they will drift from the code over time. The intended fix is to generate this
> table from the Snakefile in CI so it cannot go stale.
> *(maintainer to-do: add that generation step; until then, treat
> `workflow/Snakefile` as authoritative if the two disagree.)*

## How to read a target path

A target path is a sequence of directory and file names, each of which is either
a fixed word or a wildcard filled in from what you request. For example:

```
results/wales-latest_filter-road-primary/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
```

reads, wildcard by wildcard, as:

| Segment                     | Wildcard              | Meaning                                       |
| --------------------------- | --------------------- | --------------------------------------------- |
| `wales-latest`              | `DATASET`             | the OpenStreetMap extract                     |
| `filter-road-primary`       | `FILTER_SLUG`         | the network filter (primary roads)            |
| `hazard-aqueduct-river`     | `HAZARD_SLUG`         | the hazard dataset                            |
| `EAD_and_cost_per_RP`       | `DIRECT_DAMAGE_TYPES` | expected annual damages and cost per RP       |
| `agg-sum`                   | `AGG_FUNC_SLUG`       | aggregation function (sum)                    |
| `admin-level-0`             | `ADMIN_SLUG`          | aggregate to country level                    |

Because open-gira works backwards from the file you name, a small change in the
path — `admin-level-0` to `admin-level-1`, `wales-latest` to `kenya-latest` —
asks a materially different question, and snakemake will plan the work
accordingly. Preview it with `snakemake -n` before running.

## The wildcards

The following constraints govern what each wildcard may contain. The "pattern"
column is the regular expression from `workflow/Snakefile`; the "meaning" column
is the plain-language version.

### Places, networks, hazards

| Wildcard        | Pattern                        | Meaning                                                                 |
| --------------- | ------------------------------ | ----------------------------------------------------------------------- |
| `DATASET`       | `[^_/]+`                       | Source dataset name (no `_` or `/`), e.g. `wales-latest`.               |
| `FILTER_SLUG`   | `filter-[^_/]+`                | Network filter, always prefixed `filter-`, e.g. `filter-road-primary`. |
| `HAZARD_SLUG`   | `hazard-[^_/]+` or `nbs-[^_/]+`| A hazard dataset (`hazard-…`) or nature-based-solution layer (`nbs-…`). |
| `PROJECT_SLUG`  | `project-[^_/]+`               | A named project grouping, prefixed `project-`.                          |

### Slicing and chunking

| Wildcard       | Pattern            | Meaning                                                              |
| -------------- | ------------------ | ------------------------------------------------------------------- |
| `SLICE_SLUG`   | `slice-[0-9]+`     | One spatial slice of a sliced dataset, e.g. `slice-0`.              |
| `CHUNK_SLUG`   | `chunk-[\d]+`      | One chunk of a chunked computation, e.g. `chunk-3`.                 |
| `SAMPLE`       | `\d+`              | A sample index within a storm track set.                            |

### Direct damages and aggregation

| Wildcard              | Pattern                                                                  | Meaning                                                   |
| --------------------- | ------------------------------------------------------------------------ | -------------------------------------------------------- |
| `DIRECT_DAMAGE_TYPES` | `fraction_per_RP` \| `cost_per_RP` \| `EAD` \| `EAD_and_cost_per_RP` \| `EAD_and_cost_per_trigger` | Which damage product to compute. |
| `COST_OR_FRACTION`    | `cost` \| `fraction`                                                     | Report monetary cost, or damage fraction.                |
| `AGG_FUNC_SLUG`       | `agg-sum`                                                                | Aggregation function (currently sum).                    |
| `ADMIN_SLUG`          | `admin-level-[0-4]`                                                      | Administrative level to aggregate to (0 = country).      |

### Tropical cyclones

Tropical-cyclone analysis draws on several track sets, each with its own naming.

| Wildcard             | Pattern                                                                  | Meaning                                                        |
| -------------------- | ------------------------------------------------------------------------ | ------------------------------------------------------------- |
| `STORM_SET`          | `(?:IBTrACS\|STORM\|IRIS\|CHAZ\|emanuel)[^/]*`                            | The storm track set (historic or synthetic).                  |
| `STORM_BASIN`        | `EP\|NA\|NI\|SI\|SP\|WP`                                                  | Ocean basin (e.g. `NA` = North Atlantic).                     |
| `STORM_RP`           | `[0-9]+`                                                                  | Return period in years.                                       |
| `STORM_MODEL`        | `constant\|CMCC-CM2-VHR4\|CNRM-CM6-1-HR\|EC-Earth3P-HR\|HadGEM3-GC31-HM` | Climate model underlying a synthetic set.                     |
| `CHAZ_MODEL`         | `GCM-CESM2\|GCM-CNRM-CM6-1\|GCM-EC-Earth3\|GCM-IPSL-CM6A-LR\|GCM-MIROC6\|GCM-UKESM1-0-LL` | GCM for the CHAZ track set.                   |
| `CHAZ_SCENARIO`      | `CHAZ_SSP-[0-9]+_GCM-…_epoch-[0-9]+`                                      | A CHAZ scenario: SSP × GCM × epoch.                           |
| `EMANUEL_SCENARIO`   | `emanuel_ssp-[0-9]+_gcm-…_epoch-[0-9]+`                                   | An Emanuel-track scenario: SSP × GCM × epoch.                 |
| `IRIS_SCENARIO`      | `PRESENT\|SSP1-2050\|SSP2-2050\|SSP5-2050`                                | An IRIS scenario.                                             |
| `EVENTS_OR_RASTERS`  | `events\|wind_speed_raster\|RP_raster`                                    | Which cyclone product form to produce.                        |

### Rasters and files

| Wildcard      | Pattern                          | Meaning                                                       |
| ------------- | -------------------------------- | ------------------------------------------------------------ |
| `FILENAME`    | `[^/]+`                          | A file name (no `/`).                                         |
| `TIFF_FILE`   | `[^/.\s]+\.[tT][iI][fF][fF]?`    | A GeoTIFF file (`.tif` or `.tiff`, any case).                |
| `OUTPUT_DIR`  | `^.*results`                     | The output directory — **must end in `results`** (see below).|

## The `OUTPUT_DIR` rule

`OUTPUT_DIR` is constrained to end in the literal string `results`. This is
deliberate: it stops snakemake from matching *past* the results directory into
unrelated folders when it resolves a path, while still letting you choose where
results live (`./results`, `/data/open-gira/results`, `/my-results`, and so on
are all valid). The practical consequence — that your output directory name must
end in `results` — is recorded in
[constraints and gotchas](constraints-and-gotchas.md).

## Worked examples

A few complete targets, to show the grammar in use:

```bash
# The primary road network of Wales
results/wales-latest_filter-road-primary/edges.gpq

# Rasterised (hazard-split) road network for Egypt, one slice
results/splits/egypt-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet

# Expected annual damages, per admin-0 region, for Egyptian roads and river flooding
results/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet

# Maximum wind field for Puerto Rico's grid region, from IBTrACS storms
results/power/by_country/PRI/storms/IBTrACS/max_wind_field.nc

# Grid exposure to a single storm (Irma, 2017) across affected countries
results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/exposure_by_target.nc
```

For the workflow-by-workflow catalogue of output paths, see the
[usage guide](../user-guide/usage.md); for the rules that produce them, read
`workflow/Snakefile` and the `.smk` files it includes.
