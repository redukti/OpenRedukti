# Create base Linux package with all dependencies we need 
# for the OpenRedukti server
FROM ubuntu:18.04

LABEL description="redukti-ubuntu-base" maintainer="redukti@gmail.com"

RUN set -x \
    && apt-get update \
    && apt-get install -y libreadline-dev libopenblas-dev libsnappy-dev libbz2-dev \
    && apt-get install -y libuv1-dev librocksdb-dev libprotobuf-dev protobuf-compiler \
    && apt-get install -y build-essential autoconf libtool pkg-config \
    && apt-get install -y libgflags-dev libgtest-dev git curl wget \
    && apt-get clean

## Sadly we have to build grpc from source
RUN set -x \
    && mkdir /sources \
    && cd /sources \
    && git clone -b $(curl -L https://grpc.io/release) https://github.com/grpc/grpc \
    && cd grpc \
    && git submodule update --init \
    && cd /sources/grpc \
    && prefix=/Software/grpc make install \
    && rm -rf /sources

ENV PATH=/Software/grpc/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/grpc/lib:${LD_LIBRARY_PATH}