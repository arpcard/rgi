"""
from this directory run:

pip install .

or

python setup.py buil1d
python setup.py test
python setup.py install

"""

from app.settings import SOFTWARE_VERSION
from setuptools import setup

with open('requirements.txt') as fh:
    requirements = fh.read().splitlines()

files = ["*.py", "*.txt", "*.md", "*.json", "_data/*", "_db/*", "diff/*" ]

setup(name = "RGI",
    version = SOFTWARE_VERSION,
    description = "Resistance Gene Identifier",
    author = "Alcock et al. 2020 CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic Acid Research, 48, D517-525",
    author_email = "card@mcmaster.ca",
    url = "https://card.mcmaster.ca",
    packages = ['app'],
    package_data = {'app' : files},
    include_package_data=True,
    scripts = ["rgi"],
    long_description = """This tool provides a preliminary annotation of your DNA sequence(s) based upon the data available in CARD. Hits to genes tagged with Antibiotic Resistance ontology terms will be highlighted. As CARD expands to include more pathogens, genomes, plasmids, and ontology terms this tool will grow increasingly powerful in providing first-pass detection of antibiotic resistance associated genes.""",
    license = "The license is located at https://card.mcmaster.ca/about",
    download_url = "https://card.mcmaster.ca/download",
    install_requires=requirements,
    tests_require=requirements,
    platforms = "Linux, Mac OS X"
)
