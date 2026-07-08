# trackset

Tropical cyclone track sets, normalised.

`trackset` reads observed and synthetic tropical cyclone track sets into a
single, validated [GeoDataFrame](https://geopandas.org) schema, and carries
the metadata you need to actually use them for risk analysis — the time span
a set represents, whether its timestamps are real or fabricated, the wind
speed averaging period — alongside the data, including through GeoParquet
serialisation.

Supported sources:

| Source | Type | Distribution | Reader |
|---|---|---|---|
| [IBTrACS](https://www.ncei.noaa.gov/products/international-best-track-archive) | observed | open (CSV) | `read_ibtracs` |
| [STORM](https://doi.org/10.1038/s41597-020-0381-2) (Bloemendaal et al.) | synthetic | open, Zenodo (CSV) | `read_storm` |
| [IRIS](https://doi.org/10.1038/s41597-024-03250-y) (Sparks & Toumi) | synthetic | open (text) | `read_iris` |
| [CHAZ](https://doi.org/10.1002/2017MS001186) (Lee et al.) | synthetic | on request (netCDF) | `read_chaz_netcdf` * |
| Emanuel / WindRiskTech | synthetic | on request (MATLAB) | `read_emanuel_mat` * |

\* Raw CHAZ and Emanuel ingest is absorbed from
[chaz-track-parser](https://github.com/thomas-fred/chaz-track-parser) and
[emanuel-track-parser](https://github.com/thomas-fred/emanuel-track-parser):
these two sets need frequency calibration against observations before use
(see below). `read_chaz` / `read_emanuel` remain for GeoParquet already
produced by those external tools.

## Why

Every one of these sources uses different file formats, units, wind speed
averaging periods, missing-data conventions and identifier schemes. Every
project that consumes more than one of them rewrites the same normalisation
code — [open-gira](https://github.com/nismod/open-gira) alone carried five
separate parser scripts. `trackset` is that normalisation, extracted, tested
and importable, so that downstream work (wind field estimation, exposure,
damage and disruption analysis — in open-gira, [CLIMADA](https://climada.ethz.ch),
or your own code) can start from one schema.

## Install

```bash
pip install nismod-trackset
```

(Not yet published; for now: `pip install -e libs/trackset` from an open-gira
checkout.)

## Quickstart

```python
import shapely
import trackset

# one 1,000-year STORM sample, two basins
ts = trackset.read_storm(
    [
        "STORM_DATA_IBTRACS_NA_1000_YEARS_0.csv",
        "STORM_DATA_IBTRACS_EP_1000_YEARS_0.csv",
    ],
    sample=0,
)

ts.n_tracks           # e.g. 27_000 tracks...
ts.years              # ...representing 1,000.0 simulated years
ts.annual_frequency   # so 27.0 tracks per year
ts.data               # the GeoDataFrame of track points

# subset to tracks passing near Puerto Rico, keeping each track's points
# from first arrival to last departure
puerto_rico = shapely.box(-67.5, 17.6, -65.1, 18.6)
nearby = ts.subset(puerto_rico, buffer_deg=2.0)

# GeoParquet round-trip preserves the metadata, not just the points
nearby.to_parquet("PRI_tracks.geoparquet")
again = trackset.TrackSet.read_parquet("PRI_tracks.geoparquet")
again.source, again.years, again.synthetic_time
```

## Raw ingest and frequency calibration (CHAZ, Emanuel)

CHAZ and Emanuel tracks arrive as ragged arrays (netCDF and MATLAB) and,
unlike STORM/IRIS, are not generated at observed storm rates -- they must be
calibrated before their frequencies mean anything. The pipeline, previously
spread across two external repos, is now four calls:

```python
import numpy as np
import trackset
from trackset import frequency
from trackset.readers import filter_epoch, read_chaz_netcdf

# 1. ingest: unravel the datacube; infer radius-to-max-winds (Willoughby
#    2004) and minimum pressure (Holland 1980 pressure profile with a
#    Vickery & Wadhera 2008 shape-parameter fit); tag basins
raw = read_chaz_netcdf("CHAZ_..._sample-000.nc", genesis_method="CRH", sample=0)

# 2. window to epochs (2050 +/- 15 years, and the set's own baseline)
baseline = filter_epoch(raw, epoch=2010, half_width_years=15)
target = filter_epoch(raw, epoch=2050, half_width_years=15)

# 3. target per-basin frequency: observed rate (IBTrACS, 2002-2023) scaled
#    by the synthetic set's own relative epoch change -- trust the generator
#    for the climate-change signal, not the absolute rate
observed = frequency.observed_frequency(trackset.read_ibtracs("ibtracs.csv").data)
target_rate = observed * frequency.relative_frequency(baseline, target)

# 4. re-assign tracks to years matching that rate; the calibration tells
#    you how many years the set now represents
calibrated, years = frequency.calibrate_frequency(
    target, target_rate, rng=np.random.default_rng(0)
)
ts = frequency.finalise(
    calibrated,
    source="CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050",
    years=years,
    attributes={"ssp": 585, "gcm": "UKESM1-0-LL", "epoch": 2050},
)
```

Emanuel is the same shape, starting from
`read_emanuel_mat(path)` per input basin (no structure inference needed --
radius and pressure are provided).

Note that `calibrate_frequency` requires a *seeded* random generator: the
year re-assignment is stochastic, and the upstream implementations drew from
an unseeded global state, making their outputs unreproducible from source.
Keep the seed with your provenance metadata.

Supporting pieces, importable separately: `trackset.basins` (STORM basin
polygons and point tagging), `trackset.physics` (Willoughby 2004 RMW,
Holland 1980 pressure, Vickery & Wadhera 2008 B, Coriolis, per-basin
environmental pressures).

## The schema

A track set's `data` is a GeoDataFrame of track points with a
`DatetimeIndex`, point geometry in EPSG:4326 (longitudes in [-180, 180)),
and one row per (track, timestep). `trackset.validate` checks conformance
and reports *all* problems at once; readers validate on construction.

Required columns:

| column | type | meaning |
|---|---|---|
| `track_id` | str | unique per track within a set |
| `timestep` | int | contiguous and increasing within a track |
| `max_wind_speed_ms` | float | 1-minute sustained wind at 10m, m/s |
| `min_pressure_hpa` | float | eye pressure, hPa |
| `radius_to_max_winds_km` | float | eye to maximum winds distance, km |

Optional columns (preserved where the source provides them; consumers must
not require them): `year`, `month`, `tc_number`, `basin_id`, `sample`,
`category`, `landfall`, `distance_to_land_km`, `name`.

Readers normalise on ingest:

- **winds to 1-minute sustained**: STORM's 10-minute winds are rescaled;
  IBTrACS per-agency reports are adjusted using the Knapp & Kruk (2010)
  factors (as used by CLIMADA) before averaging across agencies
- **units**: knots to m/s, nautical miles to km
- **longitudes** to [-180, 180)
- **synthetic time**: sets that report only timestep numbers get a fabricated
  `DatetimeIndex` at the set's native frequency (flagged by
  `TrackSet.synthetic_time` — do not treat these dates as calendar time;
  they exist so tracks can be interpolated and translation speeds measured)
- **duplicate track points** (present in raw STORM/IRIS) are dropped

## Metadata

The information you need to turn a bag of tracks into frequencies travels
with the data:

- `source` — set family and variant, e.g. `"STORM-constant"`,
  `"CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050"`
- `years` — total years of observation or simulation the set represents:
  the denominator for `annual_frequency` and any expected-annual metric
  downstream. Known per sample for STORM/IRIS (1,000), inferred from the
  observed span for IBTrACS, and required from the caller for CHAZ/Emanuel
  (whose raw formats don't state it)
- `synthetic_time`, `wind_averaging_period` — see above
- `attributes` — free-form provenance (GCM, SSP, epoch, upstream filenames)

`TrackSet.to_parquet` embeds this as JSON in the parquet file's key-value
metadata; `TrackSet.read_parquet` restores it. The files remain plain
GeoParquet, readable by any tool.

## Design notes

- **Faithful first.** Readers are ports of open-gira's parsers, preserving
  its conventions (e.g. `year` offset by `1000 * sample` for STORM/IRIS so
  concatenated samples keep unique years) so open-gira can adopt the library
  without changing results. Deviations are deliberate and small: agency wind
  averaging uses pandas `mean(skipna=True)` rather than `numpy.nanmean` row
  application (equivalent, much faster); track lengths are counted by
  `groupby` rather than hashing ids; one exact knots-to-m/s constant
  (1852/3600) replaces the two truncations used upstream; the plausible
  pressure clamp lives inside `p_min_holland_1980` rather than at its call
  site; and frequency calibration takes a mandatory seeded generator.
- **Functions over frameworks.** A `TrackSet` is a thin, frozen container;
  everything real is a plain function on GeoDataFrames in `trackset.ops`
  (`saffir_simpson_category`, `wrap_longitude`, `subset_by_geometry`, ...).
  Use the container, or don't.
- **Validation is loud and complete**: `SchemaError` lists every problem,
  not just the first.

## Roadmap

- `trackset.interop.to_climada` — export to `climada.hazard.TCTracks`
  (stub documents the intended mapping)
- validate absorbed CHAZ/Emanuel ingest against the published calibrated
  track sets (requires the on-request raw data; unit tests currently cover
  fabricated miniatures of each array layout)
- slim chaz-track-parser / emanuel-track-parser workflows to thin calls
  into this package
- track interpolation to arbitrary frequency (currently lives in open-gira's
  wind field estimation; belongs here)
- adopt in open-gira: replace `workflow/tropical-cyclone/parse_*.py` with
  calls into this package

## Development

```bash
pip install -e ".[test]"
pytest
ruff format --check . && ruff check .
```

Part of [open-gira](https://github.com/nismod/open-gira) (MIT licence).
Developed at the University of Oxford. The track data itself is subject to
each provider's own licence and access terms — this package ships no data.
