#!/usr/bin/env bash
# run unit tests for RGI

# exit on failure of any command
set -e

# check if rgi is installed
rgi --help
diamond --version
bowtie2 --version
bedtools --help
samtools --help
bamtools --help
jellyfish --help
prodigal --help
# bwa

# get latest card database
wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data
mkdir -p card_data
tar xvf card_data.tar.bz2 -C card_data

# get latest card variants
wget -O prevalence-v3.0.4.tar.gz https://card.mcmaster.ca/download/6/prevalence-v3.0.4.tar.gz
mkdir -p card_variants 
tar xvf prevalence-v3.0.4.tar.gz -C card_variants

# create fasta files with annotations from card.json
rgi card_annotation --input card_data/card.json

data_version=`echo card_database_v*.fasta | sed 's/.*card_database_v\(.*\).fasta/\1/'`
variants_version=`echo prevalence-v*.tar.gz | sed 's/.*prevalence-v\(.*\).tar.gz/\1/'`

# create fasta files with annotations from variants
rgi wildcard_annotation --input_directory card_variants --version "$variants_version" --card_json card_data/card.json

# clean
rgi clean --debug

# load
rgi load --card_json card_data/card.json --card_annotation card_database_v*.fasta --wildcard_index card_variants/index-for-model-sequences.txt --wildcard_version "$variants_version" --wildcard_annotation wildcard_database_v*.fasta --debug

# check database
rgi database -v --all

# for test_1.py
cp card_data/card.json app/_data

# for test_3.py
cp card_data/card.json tests/inputs

# run unit tests
cd tests
pytest -v -rxs

# exit with the exitcode thrown by pytest
exit $?
