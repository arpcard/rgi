#!/usr/bin/env bash
# run unit tests for RGI

# exit on failure of any command
set -e

cmd="COMMAND"

echo "=================================== CHECK DEPENDENCIES VERSIONS ==================================="

blastp -version
bowtie2 --version
diamond --version
samtools --version
bamtools --version
bedtools --version
# bwa
jellyfish --version
kma -v
prodigal -v

echo "=================================== RGI EXECUTABLE LOCATION ==================================="
which rgi
rgi main --version
which python
which python3

# echo "=================================== RGI AUTO LOAD ==================================="
# rgi auto_load --clean

echo "=================================== DOWNLOAD CARD CANONICAL DATA ==================================="
# get latest card database
wget --quiet -O card_data.tar.bz2 --no-check-certificate https://card.mcmaster.ca/latest/data
mkdir -p card_data
tar xf card_data.tar.bz2 -C card_data

echo "=================================== DOWNLOAD CARD VARIANTS DATA ==================================="
# get latest card variants
wget --quiet -O prevalence-v4.0.0.tar.bz2 --no-check-certificate https://card.mcmaster.ca/download/6/prevalence-v4.0.0.tar.bz2
mkdir -p card_variants
tar xf prevalence-v4.0.0.tar.bz2 -C card_variants
gunzip card_variants/*.gz

echo "=================================== CARD CANONICAL ANNOTATIONS ==================================="
# create fasta files with annotations from card.json
echo "$cmd rgi card_annotation --input card_data/card.json"
python3 ./rgi card_annotation --input card_data/card.json

echo "=================================== VERSIONS ==================================="
# remove file with '_all.fasta'
rm card_database_v*_all.fasta
data_version=`echo card_database_v*.fasta | sed 's/.*card_database_v\(.*\).fasta/\1/'`
echo "$cmd data_version: $data_version"
variants_version=`echo prevalence-v*.tar.bz2 | sed 's/.*prevalence-v\(.*\).tar.bz2/\1/'`
echo "$cmd variants_version: $variants_version"

echo "=================================== CARD VARIANTS ANNOTATIONS ==================================="
# create fasta files with annotations from variants
echo "$cmd rgi wildcard_annotation --input_directory card_variants --version '$variants_version' --card_json card_data/card.json"
python3 ./rgi wildcard_annotation --input_directory card_variants --version "$variants_version" --card_json card_data/card.json

echo "=================================== CLEAN OLD DATABASES ==================================="
# clean
echo "$cmd rgi clean --debug"
python3 ./rgi clean --debug

echo "=================================== LOAD DATABASES ==================================="
# load
echo "$cmd rgi load --card_json card_data/card.json --card_annotation card_database_v${data_version}.fasta --wildcard_index card_variants/index-for-model-sequences.txt --wildcard_version '$variants_version' --wildcard_annotation wildcard_database_v${variants_version}.fasta --debug"
python3 ./rgi load --card_json card_data/card.json --card_annotation card_database_v${data_version}.fasta --wildcard_index card_variants/index-for-model-sequences.txt --wildcard_version "$variants_version" --wildcard_annotation wildcard_database_v${variants_version}.fasta --debug


echo "=================================== CHECK LOADED DATABASES ==================================="
# check database
echo "$cmd rgi database -v --all"
python3 ./rgi database -v --all

echo "=================================== COPY DATA & INPUTS ==================================="
# for test_1.py
echo "$cmd cp card_data/card.json app/_data"
cp card_data/card.json app/_data

# for test_3.py
echo "$cmd cp card_data/card.json tests/inputs"
cp card_data/card.json tests/inputs

echo "=================================== CLEAN ==================================="

# clean up files
rm -r card_*
rm wildcard_database_v*
rm prevalence-v*

echo "=================================== RUN TESTS ==================================="

pushd tests
echo `pwd`
tree -L 2
#pytest test_1.py -v -rxs --color=auto --durations=0 -k "protein"
pytest --capture=fd -v -rxs
popd

echo "=================================== DONE ==================================="

# exit with the exitcode thrown by pytest
exit $?
