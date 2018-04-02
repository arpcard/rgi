"""
from this directory run:

pip3 install .

or

python3: setup.py buil1d
python3: setup.py test
python3: setup.py install

"""

from app.settings import SOFTWARE_VERSION
from setuptools import setup

files = ["*.py", "*.txt", "*.md", "*.json", "_data/*", "_db/*", "diff/*" ]

setup(name = "RGI",
    version = SOFTWARE_VERSION,
    description = "Resistance Gene Identifier",
    author = "Jia et al. 2017. CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 45, D566-573.",
    author_email = "card@mcmaster.ca",
    url = "https://card.mcmaster.ca",
    packages = ['app'],
    package_data = {'app' : files},
    include_package_data=True,
    scripts = ["rgi"],
    long_description = """This tool provides a preliminary annotation of your DNA sequence(s) based upon the data available in CARD. Hits to genes tagged with Antibiotic Resistance ontology terms will be highlighted. As CARD expands to include more pathogens, genomes, plasmids, and ontology terms this tool will grow increasingly powerful in providing first-pass detection of antibiotic resistance associated genes.""",
    license = "The license is located at https://card.mcmaster.ca/about",
    download_url = "https://card.mcmaster.ca/download",
    install_requires=['biopython', 'filetype', 'pandas', 'pytest', 'mock'],
    platforms = "Linux, Mac OS X"
)
