FROM nvidia/cuda:9.0-devel-ubuntu16.04

RUN apt-get update && apt-get install --no-install-recommends -y \
        libopenblas-dev \
        python \
        python-numpy \
        gfortran \
        git \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN export CMAKE_VERSION=3.10.3 && \
    export CMAKE_VERSION_SHORT=3.10 && \
    export CMAKE_URL=https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_SCRIPT=cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    export CMAKE_PREFIX=/usr/local && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_SCRIPT} && \
    mkdir -p ${CMAKE_PREFIX} && \
    sh ${CMAKE_SCRIPT} --skip-license --prefix=${CMAKE_PREFIX} && \
    rm ${CMAKE_SCRIPT}
