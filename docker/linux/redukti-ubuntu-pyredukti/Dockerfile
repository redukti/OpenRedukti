### Runs Jupyter
FROM ubuntu:18.04

LABEL description="redukti-ubuntu-pyredukti" maintainer="redukti@gmail.com"

RUN set -x \
    && apt-get update \
    && apt-get install -y libreadline-dev libopenblas-dev \
    && apt-get install -y libprotobuf-dev protobuf-compiler \
    && apt-get install -y build-essential \
    && apt-get install -y git wget \
    && apt-get clean \
    && mkdir -p /Software \
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
    && /Software/cmake-3.14.5-Linux-x86_64/bin/cmake -DCMAKE_INSTALL_PREFIX=/Software/redukti -DCMAKE_BUILD_TYPE=Release .. \
    && make install \
    && rm -rf /Software/cmake-3.14.5-Linux-x86_64 \
    && rm -rf /sources

ENV PATH /Software/redukti/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/redukti/lib:${LD_LIBRARY_PATH}

RUN set -x \
    && apt-get install -y python3 python3-pip \
    && pip3 install cython \
    && pip3 install grpcio grpcio-tools \
    && pip3 install jupyter matplotlib

RUN set -x \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/PyRedukti.git \
    && cd PyRedukti \
    && sed -i 's/~//g' setup.py \
    && python3 setup.py bdist_wheel \
    && pip3 install dist/*.whl \
    && mkdir -p /testdata \
    && cp -r testdata/* /testdata \
    && rm -rf /sources

RUN set -x \
    && mkdir -p /data

WORKDIR /data

EXPOSE 8888
CMD jupyter notebook --allow-root --ip=0.0.0.0
