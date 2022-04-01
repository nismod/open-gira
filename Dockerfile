FROM ubuntu:20.04

# Silence interactive installations
ENV DEBIAN_FRONTEND=noninteractive

# Install required software tools
RUN apt update
RUN apt install -y \
    python3 \
    python3-dev \
    python3-pip \
    python-is-python3
RUN apt install -y \
    build-essential \
    cmake \
    gdal-bin \
    git \
    libboost-program-options-dev \
    libexpat1-dev \
    libbz2-dev \
    libgeos-dev \
    unzip \
    wget \
    zlib1g-dev

# git install remaining tools
RUN git clone https://github.com/mapbox/protozero && \
    git clone --branch v2.17.1 --depth 1 https://github.com/osmcode/libosmium && \
    git clone --branch v1.13.2 --depth 1 https://github.com/osmcode/osmium-tool

# Build osmium
RUN cd osmium-tool && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cd ../..
# Create a symlink so we can access it from anywhere on the command line
RUN cp -s /osmium-tool/build/osmium /usr/bin/osmium

# Clone the repository
COPY . open-gira

# Install Python requirements with pip
RUN pip install -r open-gira/requirements.txt

WORKDIR /open-gira

# Showtime
#RUN cd open-gira && snakemake --cores all -R all
#RUN cd open-gira && snakemake --cores all -R clean
