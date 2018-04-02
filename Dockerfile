# base image (stripped down ubuntu for Docker)
FROM phusion/baseimage

# metadata
LABEL base.image="ubuntu"
LABEL version="1"
LABEL software="RGI"
LABEL software.version="4.0.3"
LABEL description="Tool to identify resistance genes using the CARD database"
LABEL website="https://card.mcmaster.ca/"
LABEL documentation="https://github.com/arpcard/rgi/blob/master/README.rst"
LABEL license="https://github.com/arpcard/rgi/blob/master/LICENSE"
LABEL tags="Genomics"

# maintainer
MAINTAINER Finlay Maguire <finlaymaguire@gmail.com>

# install system dependencies
RUN \
    apt-get update && \
    apt-get install -y git python3 python3-dev python3-pip ncbi-blast+ prodigal wget && \
    wget http://github.com/bbuchfink/diamond/releases/download/v0.8.36/diamond-linux64.tar.gz && \
    tar xvf diamond-linux64.tar.gz && \
    mv diamond /usr/bin

# install and test rgi
RUN git clone https://github.com/arpcard/rgi
WORKDIR rgi/
RUN pip3 install -r requirements.txt && \
    pip3 install . && \
    bash test.sh

# sort the database path issue
RUN pwd
WORKDIR /usr
RUN cp /rgi/card_data/card.json /usr/local/lib/python3.5/dist-packages/app/_data/card.json

# move to workdir 
WORKDIR /data/

# set rgi executable as cmd to allow overriding
CMD ["rgi"]

