# Running on a cluster (SLURM)

Some open-gira analyses — a global run, a large storm set, many return periods —
are too large for a single machine. Since v0.4.2, open-gira supports running on a
[SLURM](https://slurm.schedmd.com/) cluster, distributing rule jobs across
compute nodes.

This page explains how that works and gives a worked example using the profile
for the University of Oxford's ARC cluster, which ships in the repository at
`config/ARC/`.

## How it works

snakemake supports cluster execution through *executor plugins*. open-gira's
environment includes
[`snakemake-executor-plugin-slurm`](https://github.com/snakemake/snakemake-executor-plugin-slurm),
so snakemake can submit each rule job to SLURM as a separate job with its own
resource request, rather than running everything in one process.

You control this with a **profile** — a directory containing a `config.yaml`
that sets the executor, the resource defaults, and the partitions available.
When you pass `--workflow-profile <dir>`, snakemake reads that configuration and
submits jobs accordingly.

The key idea is that most of open-gira works unchanged: you still request the
same target files. What changes is *where* the work runs.

## The bundled ARC profile

`config/ARC/config.yaml` is a worked example. Its important settings are:

```yaml
executor: "slurm"          # submit jobs to SLURM
jobs: 128                  # up to 128 jobs in flight at once
rerun-incomplete: true     # re-run jobs left incomplete by an interrupted run
keep-going: true           # don't stop the whole run on one failed job

default-resources:
  slurm_partition: "short"
  mem_mb: "16000"          # 16 GB per job unless a rule asks for more
  runtime: "8h"

partitions:                # the cluster's partitions and their limits
  short:
    max_runtime: 720       # minutes
    max_mem_mb: 384000
    ...
```

Two settings deserve comment:

- `default-resources` applies to any rule that does not specify its own
  resources. Individual rules can request more memory or time, and snakemake
  will honour that up to the partition limits.
- `partitions` describes the cluster to snakemake so it can choose an
  appropriate partition for a job's resource request. These values are specific
  to ARC; on another cluster they must be changed.

## Running against it

From the repository root, having loaded the environment, request a target as
usual but point snakemake at the profile:

```bash
snakemake --workflow-profile config/ARC -- <your target file>
```

snakemake will submit the necessary rule jobs to SLURM, wait for them, and
assemble the result. As always, do a dry run first:

```bash
snakemake -n --workflow-profile config/ARC -- <your target file>
```

## Adapting the profile to your cluster

The ARC profile is unlikely to match your cluster exactly. To adapt it, copy the
directory and edit:

- **`partitions`** — replace with your cluster's partition names and their
  runtime, memory, CPU, and node limits.
- **`default-resources`** — set defaults appropriate to your typical rule, in
  particular `slurm_partition`, `mem_mb`, and `runtime`.
- **`slurm-logdir`** — where SLURM job logs are written (the ARC profile uses
  `./jobs/log`).
- **`slurm_extra`** — any extra `sbatch` arguments. The ARC profile uses this
  for failure email notifications; you will want to change or remove the
  hard-coded email address.
- **`slurm-no-account`** / account settings — depending on whether your cluster
  requires an account or allocation to be named.

Consult your cluster's documentation for partition names and resource limits,
and the
[snakemake SLURM plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
for the full set of options.

## A note on shared filesystems and pixi

On many clusters the home filesystem is small or slow. pixi can install
environments in a location you choose, which is useful for putting packages on
faster or larger storage:

```bash
pixi config set detached-environments $HOME/.local/share/pixi_envs
```

See [Installation](../user-guide/installation.md) for more on pixi.
