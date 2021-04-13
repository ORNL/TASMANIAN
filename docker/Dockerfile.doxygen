FROM ubuntu:20.04

RUN apt update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        build-essential \
        cmake \
        doxygen \
        graphviz \
        git \
        texlive \
        ghostscript \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
