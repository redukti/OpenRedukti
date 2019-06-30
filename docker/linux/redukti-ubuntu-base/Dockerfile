# Create base Linux package with all dependencies we need 
# for the OpenRedukti server
FROM ubuntu:18.04

LABEL description="redukti-ubuntu-base" maintainer="redukti@gmail.com"

RUN set -x \
    && apt-get update \
    && apt-get install -y libreadline-dev libopenblas-dev libsnappy-dev libbz2-dev \
    && apt-get install -y libuv1-dev librocksdb-dev libprotobuf-dev protobuf-compiler \
    && apt-get install -y build-essential autoconf libtool pkg-config \
    && apt-get install -y libgflags-dev libgtest-dev git curl wget
##     && apt-get install -y clang libc++-dev git curl

## Sadly we have to build grpc from source
RUN set -x \
    && mkdir /sources \
    && cd /sources \
    && git clone -b $(curl -L https://grpc.io/release) https://github.com/grpc/grpc \
    && cd grpc \
    && git submodule update --init

## Build and install grpc in /Software/grpc
RUN set -x \
    && cd /sources/grpc \
    && prefix=/Software/grpc make install \
    && rm -rf /sources/grpc

## We also need the latest CMake because
## the version supplied with Ubunto Bionic does not support ROOT based package finding
RUN set -x \
    && apt-get install -y wget \
    && wget -O "/Software/cmake-3.14.5-Linux-x86_64.tar.gz" "https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz" \
    && cd /Software \
    && tar xvf "cmake-3.14.5-Linux-x86_64.tar.gz" \
    && rm -rf "/Software/cmake-3.14.5-Linux-x86_64.tar.gz"

ENV PATH=/Software/cmake-3.14.5-Linux-x86_64/bin:/Software/grpc/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/grpc/lib:${LD_LIBRARY_PATH}