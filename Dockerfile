FROM ubuntu:22.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    make \
    libbamtools-dev \
    libglpk-dev \
    libboost-dev \
    zlib1g-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/squid
COPY . .

# Compile directly so container builds do not depend on host-specific Makefile paths.
RUN mkdir -p bin && \
    g++ -std=c++11 -O2 \
      -I/usr/include/bamtools \
      src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp \
      -o bin/squid \
      -lbamtools -lglpk -lz -lm

FROM ubuntu:22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libbamtools-dev \
    libglpk-dev \
    zlib1g \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY --from=build /opt/squid/bin/squid /usr/local/bin/squid
