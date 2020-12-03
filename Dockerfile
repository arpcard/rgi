# base image (stripped down ubuntu for Docker)
FROM continuumio/miniconda3

# metadata
LABEL base.image="miniconda3"
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

# get latest version of the repo
RUN git clone https://github.com/arpcard/rgi
WORKDIR rgi

# install all dependencies matching bioconda package meta.yml
RUN conda env create -f conda_env.yml

# configure conda shell
SHELL ["conda", "run", "-n", "rgi", "/bin/bash", "-c"]

# install RGI in the repo itself
RUN pip install .

# install databases
RUN rgi auto_load 

# Move to workdir
WORKDIR /data

# set rgi executable as cmd to allow overriding
ENTRYPOINT ["conda", "run", "-n", "rgi"]
