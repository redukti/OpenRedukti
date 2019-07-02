### Runs Jupyter
FROM redukti-ubuntu-risk-base:latest

LABEL description="redukti-ubuntu-pyredukti" maintainer="redukti@gmail.com"

RUN set -x \
    && apt-get install -y python3 python3-pip \
    && pip3 install cython \
    && pip3 install grpcio grpcio-tools \
    && pip3 install jupyter

RUN set -x \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/PyRedukti.git \
    && cd PyRedukti \
    && sed -i 's/~//g' setup.py \
    && python3 setup.py bdist_wheel \
    && pip3 install dist/*.whl \
    && rm -rf /sources

RUN set -x \
    && mkdir -p /data

WORKDIR /data

EXPOSE 8888
CMD jupyter notebook --allow-root --ip=0.0.0.0
