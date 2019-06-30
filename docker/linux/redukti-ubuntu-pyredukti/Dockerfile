### Runs Jupyter
FROM redukti-ubuntu-risk-base:latest

LABEL description="redukti-ubuntu-pyredukti" maintainer="redukti@gmail.com"

EXPOSE 8888

RUN set -x \
    && apt-get install -y python3 python3-pip \
    && mkdir -p /sources \
    && cd /sources \
    && git clone https://github.com/redukti/PyRedukti.git \
    && cd PyRedukti \
    && pip3 install cython \
    && pip3 install grpcio grpcio-tools \
    && sed -i 's/~//g' setup.py \
    && python3 setup.py bdist_wheel \
    && pip3 install dist/*.whl \
    && pip3 install jupyter

RUN set -x \
    && rm -rf /sources \
    && mkdir -p /data

WORKDIR /data

CMD jupyter notebook --allow-root --ip=0.0.0.0
