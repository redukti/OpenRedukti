### Create image with OpenRedukti binaries
### Requires various dependencies from redukti-ubuntu-base
FROM redukti-ubuntu-base:latest

LABEL description="redukti-ubuntu-risk-base" maintainer="redukti@gmail.com"

RUN set -x \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/OpenRedukti.git

RUN set -x \
    && cd /sources/OpenRedukti \
    && mkdir cppbuild \
    && cd cppbuild \
    && cmake -DGRPC_SERVER=ON -DCMAKE_INSTALL_PREFIX=/Software/redukti -DGRPC_ROOT=/Software/grpc -DCMAKE_BUILD_TYPE=Release .. \
    && make install

WORKDIR /Software/redukti/bin

ENV PATH /Software/redukti/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/redukti/lib:${LD_LIBRARY_PATH}