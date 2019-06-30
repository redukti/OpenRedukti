FROM alpine:latest

RUN set -x \
    apk add --no-cache libstdc++ && \
    apk add --no-cache --virtual .build-deps git build-base autoconf automake libtool curl make cmake

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
    && sed -i 's/\-Werror//g' third_party/boringssl/CMakeLists.txt \
    && sed -i 's/\-Werror//g' Makefile \
    && prefix=/Software/grpc make install \
    && rm -rf /sources/grpc

ENV protobuf_BUILD_TESTS=OFF
RUN set -x \
    && cd /sources \
    && wget "https://github.com/protocolbuffers/protobuf/releases/download/v3.8.0/protobuf-cpp-3.8.0.tar.gz" -O "protobuf-cpp-3.8.0.tar.gz" \
    && tar xvf "protobuf-cpp-3.8.0.tar.gz" \ 
    && cd /sources/protobuf-3.8.0 \
    && ./autogen.sh \ 
    && ./configure --prefix=/Software/protobuf \
    && make install \
    && rm -rf /sources/protobuf-3.8.0

ENV PATH=/Software/grpc/bin:/Software/protobuf/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/grpc/lib:/Software/protobuf/lib:${LD_LIBRARY_PATH}

RUN set -x \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/OpenRedukti.git

RUN set -x \
    && apk add --no-cache openblas-dev readline-dev \
    && cd /sources/OpenRedukti \
    && mkdir cppbuild \
    && cd cppbuild \
    && cmake -DGRPC_SERVER=ON -DProtobuf_ROOT=/Software/protobuf -DCMAKE_INSTALL_PREFIX=/Software/redukti -DGRPC_ROOT=/Software/grpc .. \
    && make install \
    && rm -rf /source/*

ENV PATH /Software/redukti/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/redukti/lib:${LD_LIBRARY_PATH}

RUN apk del .build-deps

WORKDIR /Software/redukti/bin
