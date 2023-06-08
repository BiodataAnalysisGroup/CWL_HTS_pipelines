FROM ubuntu:latest
# DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=noninteractive
# install dependencies
ARG PYTHON_VERSION=2.7.5
# PIP - openssl version > 1.1 may be an issue (try older ubuntu images)
RUN apt-get update \
    && apt-get install -y wget gcc make openssl libffi-dev libgdbm-dev libsqlite3-dev libssl-dev zlib1g-dev \
    && apt-get clean
WORKDIR /tmp/
# Build Python from source
RUN wget https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz \
    && tar --extract -f Python-$PYTHON_VERSION.tgz \
    && cd ./Python-$PYTHON_VERSION/ \
    && ./configure --enable-optimizations --prefix=/usr/local \
    && make && make install \
    && cd ../ \
    && rm -r ./Python-$PYTHON_VERSION*
# install R-base, git and samtools
RUN apt-get install -y --no-install-recommends build-essential r-base git samtools
# copy ROSE_utils.py and ROSE_utils.pyc
COPY ROSE_utils.py /usr/local/lib/python2.7
COPY ROSE_utils.pyc /usr/local/lib/python2.7
