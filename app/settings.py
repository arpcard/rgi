import os
import sys
import logging
import json

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO

# ====================================================================================
# FUNTIONS
# ====================================================================================

def determine_path():
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except:
        print("I'm sorry, but something is wrong.")
        print("There is no __file__ variable. Please contact the author.")
        sys.exit()

# ====================================================================================
# FILEPATHS
# ====================================================================================

script_path = determine_path()

path = os.getenv('DB_PATH', os.path.join(script_path, "_db/"))
data_path = os.getenv('DATA_PATH', os.path.join(script_path, "_data/"))

# ====================================================================================
# LOGGING CONFIG
# ====================================================================================
level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
# formatter = logging.Formatter('%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')
# basic log
formatter = logging.Formatter('%(levelname)s %(asctime)s : %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

LOCAL_DATABASE = os.path.join(os.getcwd(), "localDB")

APP_NAME="Resistance Gene Identifier"
SOFTWARE_VERSION = "6.0.1"
SOFTWARE_SUMMARY = 'Use the Resistance Gene Identifier to predict resistome(s) from protein or nucleotide \
data based on homology and SNP models. Check https://card.mcmaster.ca/download for software and data updates. \
Receive email notification of monthly CARD updates via the CARD Mailing List \
(https://mailman.mcmaster.ca/mailman/listinfo/card-l)'

GALAXY_PROJECT_WRAPPER='GALAXY_DATABASE must contain the following: \
data files (card.json, proteindb.fsa),\
diamond blast database (protein.db.dmnd), \
ncbi blast database (protein.db.phr, protein.db.pin, protein.db.psq)'
