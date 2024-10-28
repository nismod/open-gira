# Utilities

`open-gira` comes with a few small utilities for data processing and file management.

## Removing intermediate files

You can remove intermediate files by running the `clean` rule.

```bash
snakemake --cores 1 -R clean
```

## Geoparquet -> Geopackage

As standard we use the [GeoParquet](https://geoparquet.org/) (`.geoparquet` or
`.gpq`) format to store vector data on disk. Unfortunately common GIS software
such as QGIS may not yet support this file format. To convert file(s) to
geopackage, use:

```
python workflow/scripts/pq_to_gpkg.py <path_to_geoparquet_1> <path_to_geoparquet_2> <...>
```

This will write `.gpkg` files beside their source `.geoparquet`.

## Unpickling interactive plots

`matplotlib` plots can be interactive (zoom, pan, etc.), but not as static
images. Some rules produce pickled plot files. To view these, use:

```
python workflow/scripts/unpickle_plot.py <path_to_pickled_plot>
```

## Archiving results

The bash script `archive_results.sh` can be used to back up analysis results.
Here's an example usage:

```bash
./archive_results.sh results/ /mnt/backup/open-gira
```

This will create a timestamped folder in path of the second argument, e.g.
`/mnt/backup/open-gira/2023-07-24T101756+0100`, containing an archive of results
(excluding input files downloaded from the internet) and a README describing the
state of the repository at the time of archive creation.
