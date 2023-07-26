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
    r-base \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    gfortran \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c bioconda -c conda-forge f5c python=3.8 minimap2 samtools
RUN /opt/conda/bin/python -m pip install --upgrade pip
RUN /opt/conda/bin/python -m pip install Cython numpy
RUN /opt/conda/bin/python -m pip install m6anet
ENV PATH=/opt/conda/bin/:$PATH

RUN git clone https://github.com/hasindu2008/f5c.git
RUN /f5c/scripts/install-vbz.sh
ENV HDF5_PLUGIN_PATH=/root/.local/hdf5/lib/plugin

RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('IRanges')"  
RUN R -e "BiocManager::install('GenomicRanges')" 
RUN R -e "BiocManager::install('ensembldb')"
