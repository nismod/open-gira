FROM osgeo/proj

# Silence interactive installations
ENV DEBIAN_FRONTEND=noninteractive

# Install required software tools
RUN apt update && apt install \
    build-essential \
    cmake \
    libboost-program-options-dev \
    libexpat1-dev \
    zlib1g-dev \
    libbz2-dev \
    gdal-bin

# git install remaining tools
RUN git clone https://github.com/mapbox/protozero && \
    git clone --branch v2.18.0 --depth 1 https://github.com/osmcode/libosmium && \
    git clone --branch v1.14.0 --depth 1 https://github.com/osmcode/osmium-tool

# Build osmium
RUN cd osmium-tool && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cd ../..

# Create a symlink so we can access it from anywhere on the command line
RUN cp -s /osmium-tool/build/osmium /usr/bin/osmium

# Copy in the repository
COPY . open-gira
WORKDIR /open-gira

# Install Python requirements
RUN curl micro.mamba.pm/install.sh | bash
RUN micromamba create -f environment.yml
