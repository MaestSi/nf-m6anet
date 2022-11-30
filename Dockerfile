FROM continuumio/miniconda3

MAINTAINER Simone Maestri <simone.maestri@iit.it>

RUN apt-get update  && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    wget \
    gfortran \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c bioconda nanopolish minimap2 samtools python==3.8

RUN git clone -b read_id https://github.com/GoekeLab/m6anet.git \
&& cd m6anet \
&& python setup.py install
