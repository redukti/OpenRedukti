FROM quay.io/pypa/manylinux2010_x86_64

RUN set -x \
    && yum install -y openblas-devel wget readline-devel

RUN set -x \
    && cd / \
    && mkdir sources \
    && mkdir Software

RUN set -x \
    && wget -O "/Software/cmake-3.14.5-Linux-x86_64.tar.gz" "https://github.com/Kitware/CMake/releases/download/v3.14.5/cmake-3.14.5-Linux-x86_64.tar.gz" \
    && cd /Software \
    && tar xvf "cmake-3.14.5-Linux-x86_64.tar.gz" \
    && rm "cmake-3.14.5-Linux-x86_64.tar.gz"

RUN set -x \
    && cd /sources \
    && wget "https://github.com/protocolbuffers/protobuf/releases/download/v3.8.0/protobuf-cpp-3.8.0.tar.gz" -O "protobuf-cpp-3.8.0.tar.gz" \
    && tar xvf "protobuf-cpp-3.8.0.tar.gz" \
    && cd /sources/protobuf-3.8.0 \
    && ./autogen.sh \
    && ./configure --prefix=/Software/protobuf \
    && make install

ENV PATH /Software/protobuf/bin:${PATH}

RUN set -x \
    && cd /sources \
    && git clone -b $(curl -L https://grpc.io/release) https://github.com/grpc/grpc \
    && cd /sources/grpc \
    && git submodule update --init \
    && sed -i 's/\-Werror//g' third_party/boringssl/CMakeLists.txt \
    && sed -i 's/\-Werror//g' Makefile \
    && prefix=/Software/grpc make install

ENV PATH /Software/grpc/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/grpc/lib:/Software/protobuf/lib:${LD_LIBRARY_PATH}

ENV PATH /Software/cmake-3.14.5-Linux-x86_64/bin:${PATH}

RUN set -x \
    && cd /sources \
    && git clone https://github.com/redukti/OpenRedukti.git \
    && cd /sources/OpenRedukti \
    && mkdir cppbuild \
    && cd cppbuild \
    && cmake -DGRPC_SERVER=ON -DProtobuf_ROOT=/Software/protobuf -DCMAKE_INSTALL_PREFIX=/Software/redukti -DGRPC_ROOT=/Software/grpc .. \
    && make install

ENV PATH /Software/redukti/bin:${PATH}
ENV LD_LIBRARY_PATH /Software/redukti/lib:${LD_LIBRARY_PATH}

RUN set -x \
    && cd /sources \
    && git clone https://github.com/redukti/PyRedukti.git \
    && /opt/python/cp37-cp37m/bin/pip install cython \
    && /opt/python/cp37-cp37m/bin/pip install grpcio grpcio-tools \
    && cd /sources/PyRedukti \
    && sed -i 's/~//g' setup.py \
    && /opt/python/cp37-cp37m/bin/python3 setup.py bdist_wheel \
    && /opt/python/cp37-cp37m/bin/pip install dist/*.whl \
    && /opt/python/cp37-cp37m/bin/pip install jupyter

ENV PATH /opt/python/cp37-cp37m/bin:${PATH}

RUN set -x \
    && rm -rf /sources

RUN set -x \
    && rm -rf /Software/cmake-3.14.5-Linux-x86_64