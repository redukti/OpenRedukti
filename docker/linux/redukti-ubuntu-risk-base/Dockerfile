### Create image with OpenRedukti binaries
### Requires various dependencies from redukti-ubuntu-base
FROM redukti-ubuntu-base:latest

LABEL description="redukti-ubuntu-risk-base" maintainer="redukti@gmail.com"

ENV PATH=:/Software/grpc/bin:${PATH}

RUN set -x \
    && wget -O "/Software/cmake-3.14.5-Linux-x86_64.tar.gz" "https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz" \
    && cd /Software \
    && tar xvf "cmake-3.14.5-Linux-x86_64.tar.gz" \
    && rm -rf "/Software/cmake-3.14.5-Linux-x86_64.tar.gz" \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/OpenRedukti.git \
    && cd /sources/OpenRedukti \
    && mkdir cppbuild \
    && cd cppbuild \
    && /Software/cmake-3.14.5-Linux-x86_64/bin/cmake -DGRPC_SERVER=ON -DCMAKE_INSTALL_PREFIX=/Software/redukti -DGRPC_ROOT=/Software/grpc -DCMAKE_BUILD_TYPE=Release .. \
    && make install \
    && rm -rf /Software/cmake-3.14.5-Linux-x86_64 \
    && rm -rf /sources

WORKDIR /Software/redukti/bin

ENV PATH /Software/redukti/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/redukti/lib:${LD_LIBRARY_PATH}