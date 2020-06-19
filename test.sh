#!/usr/bin/env bash
# run unit tests for RGI

# exit on failure of any command
set -e

cmd="COMMAND"

echo "=================================== RGI EXECUTABLE LOCATION ==================================="
which rgi
rgi -h

echo "=================================== DOWNLOAD CARD CANONICAL DATA ==================================="
# get latest card database
wget -O card_data.tar.bz2 --no-check-certificate https://card.mcmaster.ca/latest/data
mkdir -p card_data
tar xf card_data.tar.bz2 -C card_data

echo "=================================== DOWNLOAD CARD VARIANTS DATA ==================================="
# get latest card variants
wget -O prevalence-v3.0.4.tar.gz --no-check-certificate https://card.mcmaster.ca/download/6/prevalence-v3.0.4.tar.gz
mkdir -p card_variants 
tar xf prevalence-v3.0.4.tar.gz -C card_variants

echo "=================================== CARD CANONICAL ANNOTATIONS ==================================="
# create fasta files with annotations from card.json
echo "$cmd python3 ./rgi card_annotation --input card_data/card.json"
python3 ./rgi card_annotation --input card_data/card.json

echo "=================================== VERSIONS ==================================="

data_version=`echo card_database_v*.fasta | sed 's/.*card_database_v\(.*\).fasta/\1/'`
echo "$cmd data_version: $data_version"
variants_version=`echo prevalence-v*.tar.gz | sed 's/.*prevalence-v\(.*\).tar.gz/\1/'`
echo "$cmd variants_version: $variants_version"

echo "=================================== CARD VARIANTS ANNOTATIONS ==================================="
# create fasta files with annotations from variants
echo "$cmd python3 ./rgi wildcard_annotation --input_directory card_variants --version '$variants_version' --card_json card_data/card.json"
python3 ./rgi wildcard_annotation --input_directory card_variants --version "$variants_version" --card_json card_data/card.json

echo "=================================== CLEAN OLD DATABASES ==================================="
# clean
echo "$cmd ./rgi clean --debug"
python3 ./rgi clean --debug

echo "=================================== LOAD DATABASES ==================================="
# load
echo "$cmd python3 ./rgi load --card_json card_data/card.json --card_annotation card_database_v${data_version}.fasta --wildcard_index card_variants/index-for-model-sequences.txt --wildcard_version '$variants_version' --wildcard_annotation wildcard_database_v${variants_version}.fasta --debug"
python3 ./rgi load --card_json card_data/card.json --card_annotation card_database_v${data_version}.fasta --wildcard_index card_variants/index-for-model-sequences.txt --wildcard_version "$variants_version" --wildcard_annotation wildcard_database_v${variants_version}.fasta --debug

echo "=================================== CHECK LOADED DATABASES ==================================="
# check database
echo "$cmd python3 ./rgi database -v --all"
python3 ./rgi database -v --all

echo "=================================== COPY DATA & INPUTS ==================================="
# for test_1.py
echo "$cmd card_data/card.json app/_data"
cp card_data/card.json app/_data

# for test_3.py
echo "$cmd card_data/card.json tests/inputs"
cp card_data/card.json tests/inputs

echo "=================================== RUN TESTS ==================================="
# run unit tests
echo "$cmd cd tests"
cd tests
echo "$cmd pytest --capture=fd -v -rxs"
pytest --capture=fd -v -rxs

pwd

rm -r ../card_*
rm ../prevalence-v*.gz
rm ../wildcard_database_v*

echo "=================================== DONE ==================================="

# exit with the exitcode thrown by pytest
exit $?

