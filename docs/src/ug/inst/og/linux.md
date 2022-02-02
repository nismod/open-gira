# Linux installation

If you have [conda](https://docs.conda.io/en/latest/) installed, it is probably
simplest to follow the [conda install instructions](mac.md#conda) in the MacOS section.

The workflow for installing on Debian/Ubuntu is the one used in the Dockerfile.
For installations on different brands of Linux, please follow the steps for
[MacOS installation with pip](mac.md#pip), and note that the osmium-tool
distributions are frequently out of date, and you will likely have to
build the osmium-tool from its source.

## Debian/Ubuntu

### Dependencies

The first step is to get hold of the packages we'll need to run the software:

```shell
sudo apt update && sudo apt install \
  build-essential \
  cmake \
  git \
  libboost-program-options-dev \
  libexpat1-dev \
  libbz2-dev \
  python3 \
  python3-pip \
  zlib1g-dev
```

Enter 'yes' to install when prompted, and address any other questions in the install process if they arise.

### osmium-tool

The next step is to install the osmium-tool and its library.
This needs to be built from source, because the version distributed by `apt` is out of date.

```shell
git clone https://github.com/mapbox/protozero && \
    git clone --branch v2.17.1 --depth 1 https://github.com/osmcode/libosmium && \
    git clone --branch v1.13.2 --depth 1 https://github.com/osmcode/osmium-tool
```

Now we've downloaded the source files we need to build them:

```shell
cd osmium-tool && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cd ../..
```

The last step to installing osmium is to create a symlink so we can access it from anywhere on the command line:

```shell
cp -s osmium-tool/build/osmium /usr/bin/osmium
```

## open-gira

We install open-gira by cloning the repository:

```shell
git clone https://github.com/nismod/open-gira.git
```

## Install Python requirements with pip

The final step is to make sure we have the Python libraries we need to run open-gira. 
We use pip to do this:

```shell
pip install -r open-gira/requirements.txt
```