import os, sys, json, csv, argparse, multiprocessing
from app.settings import *
from collections import OrderedDict
from Bio import SeqIO, Seq

class CARDkmers(object):
    """
    Queries sequences against CARD*kmers
    """

    def __init__(self, input, bwt, rgi, fasta, k, output, local, debug):
        # from arguments
        self.input_file = input
        self.bwt = bwt
        self.rgi = rgi
        self.fasta = fasta
        self.k = k
        self.output = output
        self.local_database = local

        # database
        self.db = path
        self.data = data_path

        if self.local_database:
            self.db = LOCAL_DATABASE
            self.data = LOCAL_DATABASE

        self.kmer_db = os.path.join(self.data, "{}mer_database.json".format(k))
        self.amr_kmers = os.path.join(self.data, "amr_{}mer.txt".format(k))

        # files
        if self.bwt:
            self.input_bam_file = input
        if self.rgi:
            self.input_json_file = input
        # if self.fasta:



        # outputs
        self.working_directory = os.path.join(os.getcwd())

        self.debug = debug
        if self.debug:
            logger.setLevel(10)

    def __repr__(self):
    	"""Returns CARDkmers class full object."""
    	return "CARDkmers({}".format(self.__dict__)

    def check_databases_exist(self):
        # Check if required files loaded
        if not os.path.exists(self.kmer_db):
            logger.error("Missing {}.".format(self.kmer_db))
            exit()
        if not os.path.exists(self.amr_kmers):
            logger.error("Missing {}.".format(self.amr_kmers))
            exit()

    # def extract_sequences_from_bam():
    #     os.system("samtools view {bam} | while read line; do awk '{print ">"$1"__"$3"__"$2"__"$5"\n"$10}'; done > {out}")
    #         )

    def run(self):
        # print args
        logger.info(json.dumps(self.__dict__, indent=2))

        logger.info("check for databases")
        self.check_databases_exist()

        # checks only one data type given
        if sum([self.bwt, self.rgi, self.fasta]) > 1:
            logger.error("Only specify one input type.")

        # if self.bwt:



        # logger.info("")
