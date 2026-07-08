# Docker

open-gira previously shipped a `Dockerfile` in the repository root. It has been
removed: it no longer built (it used a deprecated base image and a fragile
build-from-source step), and it predated the move to the [pixi](https://pixi.sh/)
lockfile, so it installed dependencies in a way the project no longer uses.

Rather than leave a broken image in place, this page describes how to run
open-gira in a container from a clean, current base — and, first, why you may
not need a container at all.

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

## A fresh Docker setup

The following is a minimal, current setup built on the official pixi image and
the project's own lockfile. It does no build-from-source: the geospatial stack,
including `osmium-tool` and GDAL, comes from the locked conda environment.

> This image is not yet part of the project's continuous integration, so build
> it locally to confirm it works in your environment before relying on it. It is
> constructed from the maintained pixi base image and the committed
> `pixi.lock`, so it should install exactly the dependencies CI uses.

Create a file named `Dockerfile` in the repository root:

```dockerfile
# syntax=docker/dockerfile:1

# The pixi image ships the pixi package manager on a Debian base.
# Pin the same pixi version the project's CI uses; check the available tags at
# https://github.com/prefix-dev/pixi/pkgs/container/pixi
FROM ghcr.io/prefix-dev/pixi:0.59.0

# A few workflow rules shell out to these system utilities.
RUN apt-get update \
    && apt-get install -y --no-install-recommends wget unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /open-gira

# Install the locked environment first so this layer is cached until the
# dependencies themselves change.
COPY pixi.toml pixi.lock ./
RUN pixi install --locked

# Copy the rest of the workflow.
COPY . .

# Run commands inside the pixi environment. With this entrypoint, whatever you
# pass after the image name runs in the environment; the default is an
# interactive shell.
ENTRYPOINT ["pixi", "run"]
CMD ["bash"]
```

Add a `.dockerignore` alongside it so local artefacts do not bloat the build
context or leak into the image:

```gitignore
.git
.pixi
results
book
docs/book
__pycache__
*.pyc
```

### Build and run

```bash
# Build the image
docker build --network host -t open-gira .

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

> **Building behind a proxy.** If your network routes outbound traffic through a
> proxy, the build's package downloads must go through it too. Pass the proxy
> settings and, if needed, a CA certificate into the build (for example with
> `--network host` and `--build-arg` values), following your platform's guidance.

## Other options

- **Dev Container.** If you use VS Code, a
  [Dev Container](https://containers.dev/) referencing the same pixi base image
  gives an in-container development environment with the editor attached. This
  is a good fit if you want to *develop* open-gira in a container rather than
  only run it.
- **Docker Compose.** The repository's `docker-compose.yml` retains a service
  for building the documentation with mdBook. Once you have created a
  `Dockerfile` from the setup above, you can re-add an `app` service that builds
  it; a commented stub is left in the compose file.
- **Your own base image.** If you already standardise on a particular
  conda/mamba or Python base image, install pixi into it and run
  `pixi install --locked` as above. The essential ingredients are the pixi
  binary and the committed `pixi.lock`.
