FROM ubuntu:18.04

RUN apt update && apt install --no-install-recommends -y \
        build-essential \
        clang \
        python \
        python-numpy \
        python3 \
        python3-numpy \
        cmake \
        gfortran \
        git \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
