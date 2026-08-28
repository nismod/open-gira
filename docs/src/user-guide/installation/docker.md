# Docker

open-gira ships a `Dockerfile` for running the workflow in a container. It is
built on the official [pixi](https://pixi.sh/) image and the project's committed
`pixi.lock`, and does no build-from-source: the geospatial stack, including
`osmium-tool` and GDAL, comes from the locked conda environment.

(An earlier Dockerfile was replaced. It no longer built — it used a deprecated
base image and a fragile build-from-source step — and predated the move to the
pixi lockfile.)

Before reaching for it, it is worth knowing you may not need a container at all.

## You may not need Docker

open-gira's supported installation path is pixi, which already solves the
problem Docker is often reached for: a reproducible environment pinned by a
lockfile. On Linux and macOS, a native install with

```bash
pixi install
pixi shell
```

gives you the same locked dependencies a container would, without the container.
See [Installation](../installation.md) and the platform notes for
[Linux / Mac](linux-mac.md) and [Windows](windows.md). If your reason for Docker
is reproducibility, pixi may already be enough.

Reach for a container when you need isolation from the host, a clean environment
on a shared machine, or a reproducible unit to run in CI or on a platform that
expects an image.

## Build and run

From a clone of the repository:

```bash
# Build the image
docker build -t open-gira .

# Open an interactive shell in the environment
docker run --rm -it open-gira

# Run a target directly (everything after the image name runs via `pixi run`)
docker run --rm -it open-gira snakemake -n --cores 2 -- \
  results/wales-latest_filter-road-primary/edges.gpq
```

To keep results on the host, mount a directory over the container's output
location:

```bash
docker run --rm -it -v "$PWD/results:/open-gira/results" open-gira \
  snakemake --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

The image is built on every push to `main` by the `Docker image` GitHub Actions
workflow (`.github/workflows/docker.yml`), so a change that breaks the build is
caught rather than discovered later.

> **Building behind a proxy.** If your network routes outbound traffic through a
> proxy, the build's package downloads must go through it too. Pass the proxy
> settings and, if needed, a CA certificate into the build (for example with
> `--network host` and appropriate build arguments), following your platform's
> guidance.

## Developing in a container

The repository also includes a [Dev Container](https://containers.dev/)
definition at `.devcontainer/devcontainer.json`. In an editor that supports it
(for example VS Code with the Dev Containers extension), "Reopen in Container"
builds the same image, mounts your working tree, installs the locked
environment, and attaches the editor — giving a ready-to-use development
environment with Python, Ruff, and Snakemake tooling configured.

This is the option to choose when you want to *develop* open-gira in a
container rather than only run it.

## Other options

- **Docker Compose.** The repository's `docker-compose.yml` defines an `app`
  service that builds this image, and an `mdbook` service for previewing the
  documentation. For example, `docker compose run --rm app` opens a shell in the
  environment.
- **Your own base image.** If you already standardise on a particular
  conda/mamba or Python base image, install pixi into it and run
  `pixi install --locked`. The essential ingredients are the pixi binary and the
  committed `pixi.lock`.
