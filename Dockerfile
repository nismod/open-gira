# A minimal, current container for open-gira, built on the official pixi image
# and the committed pixi.lock. It does no build-from-source: the geospatial
# stack, including osmium-tool and GDAL, comes from the locked conda environment.
#
# Pin the same pixi version the project's CI uses; available tags are listed at
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
# interactive shell. For example:
#   docker build -t open-gira .
#   docker run --rm -it open-gira
#   docker run --rm -it open-gira snakemake -n --cores 2 -- \
#     results/wales-latest_filter-road-primary/edges.gpq
ENTRYPOINT ["pixi", "run"]
CMD ["bash"]
