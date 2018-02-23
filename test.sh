#!/usr/bin/env bash
# run unit tests for RGI

# exit on failure of any command
set -e

# get latest card database
wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data
mkdir -p card_data
tar xvf card_data.tar.bz2 -C card_data

# load rgi
rgi load -i card_data/card.json

# for test_1.py
cp card_data/card.json app/_data

# for test_3.py
cp card_data/card.json tests/input

# run unit tests
cd tests
pytest -v -rxs

# exit with the exitcode thrown by pytest
exit $?
