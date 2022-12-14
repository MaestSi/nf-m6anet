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

RUN conda install -c bioconda f5c minimap2 samtools python==3.8
RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('IRanges')"  
RUN R -e "BiocManager::install('GenomicRanges')" 
RUN R -e "BiocManager::install('ensembldb')"

RUN git clone https://github.com/GoekeLab/m6anet.git \
&& cd m6anet \
&& python setup.py install
