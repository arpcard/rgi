# base image (stripped down ubuntu for Docker)
FROM continuumio/miniconda3

# metadata
LABEL base.image="miniuonda3"
LABEL version="2"
LABEL software="RGI"
LABEL description="Tool to identify resistance genes using the CARD database"
LABEL website="https://card.mcmaster.ca/"
LABEL documentation="https://github.com/arpcard/rgi/blob/master/README.rst"
LABEL license="https://github.com/arpcard/rgi/blob/master/LICENSE"
LABEL tags="Genomics"

# maintainer
MAINTAINER Finlay Maguire <finlaymaguire@gmail.com>

# get some system essentials
RUN apt-get update && apt-get install -y wget && conda init bash

# install rgi and system dependencies
RUN conda create --name rgi --channel conda-forge --channel bioconda rgi 

# download latest card database
RUN mkdir -p /card_data
WORKDIR /card_data
RUN wget -O data.tar.bz2 https://card.mcmaster.ca/latest/data && \
        mkdir -p canonical &&  tar xf data.tar.bz2 -C canonical

# download wildcard data (modified from @nebfield)
RUN wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants && \
    mkdir -p wildcard && \
    tar xf wildcard_data.tar.bz2 -C wildcard && \
        gunzip wildcard/*.gz

# configure conda shell
SHELL ["conda", "run", "-n", "rgi", "/bin/bash", "-c"]

# install CARD database
RUN rgi card_annotation -i canonical/card.json > card_annotation.log 2>&1 && \
    rgi load --local -i canonical/card.json --card_annotation card_database_*.fasta

# install WILDCARD database
# version number not in downloaded data so can't do this in a way that will
# automatically be entered if the wildcard download contained a version file etc
RUN rgi wildcard_annotation -i wildcard --card_json canonical/card.json -v latest > wildcard_annotation.log 2>&1 && \
    rgi load --local --wildcard_annotation wildcard_database* --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database*

# install kmer pathogen-of-origin database
RUN rgi load --local --kmer_database wildcard/61_kmer_db.json --amr_kmers wildcard/all_amr_61mers.txt --kmer_size 61 --debug > kmer_load.61.log 2>&1

# Move to workdir
WORKDIR /data

# set rgi executable as cmd to allow overriding
ENTRYPOINT ["conda", "run", "-n", "rgi"]
