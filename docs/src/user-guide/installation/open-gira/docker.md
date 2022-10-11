# Installation via Docker

It is possible to run `open-gira` within a container using Docker. You must have [Docker
installed](https://docs.docker.com/engine/install/).

1. Clone this repository

```shell
git clone https://github.com/nismod/open-gira.git
```

2. Build container image

```shell
docker build -t open-gira .
```

3. Run container

We use the `-it` options to make the container interactive and give it a terminal,
and the `-d` option to keep it running in the background (daemon mode).
We give it a nice name we can refer to later, and make sure it is created using the
image we just built.
We tell it to run `bash` so that it has something to do and doesn't close immediately.

```shell
docker run -itd --name og open-gira bash
```

4. Enter container

We can now lauch a bash shell in our new container:

```shell
docker exec -it og /bin/bash
```

And once inside we can navigate to the open-gira folder.

```shell
cd open-gira
```

We're now ready to run snakemake, launch a python prompt, or do whatever else we like inside
the container.
