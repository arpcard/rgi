"""
from this directory run:

pip3 install .

or

python3: setup.py buil1d
python3: setup.py test
python3: setup.py install

entry_points =[
	'console_scripts':['rgi = main.run']
]
"""
# from distutils.core import setup
from app.settings import SOFTWARE_VERSION
from setuptools import setup
#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.

files = ["*.py", "*.txt", "*.md", "*.json", "_data/*", "_db/*", "diff/*" ]

setup(name = "RGI",
    version = SOFTWARE_VERSION,
    description = "Resistance Gene Identifier",
    author = "Jia et al. 2017. CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 45, D566-573.",
    author_email = "card@mcmaster.ca",
    url = "https://card.mcmaster.ca",
    #Name the folder where your packages live:
    #(If you have other packages (dirs) or modules (py files) then
    #put them into the package directory - they will be found
    #recursively.)
    packages = ['app'],
    #'package' package must contain files (see list above)
    #I called the package 'package' thus cleverly confusing the whole issue...
    #This dict maps the package name =to=> directories
    #It says, package *needs* these files.
    package_data = {'app' : files},
    include_package_data=True,
    # Runners
    #'rgi' is in the root run rgi (Resistance Gene Identifier).
    #'rgi_jsontab' is in the root to convert json to tab-delimited
    #'rgi_clean' is in the root to clean blast databases
    #'rgi_load' is in the root to load new card.json file.
    scripts = ["rgi"],
    long_description = """This tool provides a preliminary annotation of your DNA sequence(s) based upon the data available in CARD. Hits to genes tagged with Antibiotic Resistance ontology terms will be highlighted. As CARD expands to include more pathogens, genomes, plasmids, and ontology terms this tool will grow increasingly powerful in providing first-pass detection of antibiotic resistance associated genes.""",
    #
    #This next part it for the Cheese Shop, look a little down the page.
    #classifiers = []
    license = "The license is located at https://card.mcmaster.ca/about",
    #location where the package may be downloaded
    download_url = "https://card.mcmaster.ca/download",
    install_requires=['biopython', 'filetype', 'pandas', 'pytest', 'mock'],
    platforms = "Linux, Mac OS X"
)
