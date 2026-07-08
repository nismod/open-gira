# Constraints and gotchas

Every real tool has sharp edges. This page collects the ones in open-gira that
have surprised users, so that they surprise you less. Where a gotcha has an open
issue, we link it; where it is enforced by the code, we say so. Some of these are
things we intend to fix rather than document forever — but until they are fixed,
knowing about them saves time.

## Your output directory must end in `results`

open-gira constrains the output directory to end in the literal string
`results` (this is the `OUTPUT_DIR` wildcard constraint in
`workflow/Snakefile`). So `./results`, `/data/open-gira/results`, and
`/scratch/my-results` are all fine, but a directory that does *not* end in
`results` will not match, and snakemake will report that it cannot produce your
target.

Why it exists: it stops snakemake, when resolving a path backwards, from
matching *past* the results directory into unrelated folders. It is a pragmatic
guard rather than a fundamental requirement. See
[Targets and wildcards](targets-and-wildcards.md) for the constraint itself.

**If you hit it:** rename your output directory so its name ends in `results`.

## The configuration file is `config.yaml`, not `config.yml`

The workflow reads `config/config.yaml`. Some older text refers to `config.yml`;
the file the code actually loads is `config.yaml`. If you copy a command from an
old note and snakemake cannot find your configuration, check the extension.

## Hazard, filter, and dataset names cannot contain `_` or `/`

The names you give hazard datasets, network filters, and infrastructure datasets
in `config/config.yaml` must not contain `_` or `/`. These characters are used as
separators when open-gira composes and decomposes target paths, so a name
containing them would make paths ambiguous. The workflow checks this at start-up
and raises a clear error if a name breaks the rule.

## `slice_count` must be a square number (or 1)

Spatial slicing partitions a bounding box into a grid, so the number of slices
must be a perfect square — `1`, `4`, `9`, `16`, `64`, and so on. Setting
`slice_count` to a non-square integer is rejected at start-up. A useful
heuristic is the first square number greater than your CPU count.

## Known issues to be aware of

The following are tracked in the issue tracker. The issues themselves are the
authoritative description; they are listed here so you know to look.

- **Relative paths in the slicing rule**
  ([#142](https://github.com/nismod/open-gira/issues/142)) — the slicing step
  has been sensitive to relative versus absolute paths. If slicing behaves
  unexpectedly, check this issue for the current status and any workaround.
- **GeoPackage modification-time trap in QGIS**
  ([#116](https://github.com/nismod/open-gira/issues/116)) — GeoPackage outputs
  can interact badly with QGIS's handling of file modification times, which can
  make snakemake think an output is stale. See the issue for details.

*(maintainer to-do: as these issues are fixed, remove them here and note the
version in which the fix landed.)*

## When in doubt, dry-run

Most confusion about "why is open-gira doing that?" is answered by asking
snakemake to explain itself before it acts:

```bash
snakemake -n --cores 2 -- <your target>
```

A dry run lists exactly which rules would fire and why, without changing
anything. It is the cheapest debugging tool open-gira offers.
