# base image (stripped down ubuntu for Docker)
FROM continuumio/miniconda3

# metadata
LABEL base.image="miniconda3"
LABEL version="1"
LABEL software="RGI"
LABEL software.version="5.1.0"
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

# install CARD database
SHELL ["conda", "run", "-n", "rgi", "card_annotation", "-i",  "canonical/card.json", ">", "card_annotation.log", "2>&1"]
SHELL ["conda", "run", "-n", "rgi", "rgi", "load", "-i", "caonical/card.json", "--card_annotation", "card_database_*.fasta"]

# install WILDCARD database
# version number not in downloaded data so can't do this in a way that will
# automatically be entered
SHELL ["conda", "run", "-n", "rgi", "rgi", "wildcard_annotation", "-i", "wildcard", "--card_json", "canonical/card.json", "-v", "docker_version_number", ">", "wildcard_annotation.log", "2>&1"]
SHELL ["conda", "run", "-n", "rgi", "rgi", "load", "--wildcard_annotation", "wildcard_database*", "--wildcard_index", "wildcard/index-for-model-sequences.txt", "--card_annotation", "card_database*"]

# install kmer pathogen-of-origin database
SHELL ["conda", "run", "-n", "rgi", "rgi", "load", "--kmer_database", "wildcard/61_kmer_db.json", "--amr_kmers", "wildcard/all_amr_61mers.txt", "--kmer_size", "61", "--debug", ">", "kmer_load.61.log", "2>&1"]

# Remove non-loaded databases to reduce container size
RUN cp *.log /data && rm -rf /card_data

WORKDIR /data
# set rgi executable as cmd to allow overriding
ENTRYPOINT ["conda", "run", "-n", "rgi", "rgi"]
