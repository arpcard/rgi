import os, sys, json, csv, argparse
from rgi.settings import *

def main(args):
    working_directory = os.getcwd()
    print(args)
    print("""
    TODO:
     - Please note that the index file used is 'baits-probes-with-sequence-info.txt'.
     - The index file should be modelled from this file (baits-probes-with-sequence-info.txt)
     - Running this command should yield a similar FASTA file like 'bait-80-20-q875t875id99.fas'
     - For now don't use this command to annotate baits, just load files
      'bait-80-20-q875t875id99.fas' and 'baits-probes-with-sequence-info.txt' using the rgi load commnand
    """)

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi baits_annotation",description='Creates baits annotations for RGI BWT from baits')
    parser.add_argument('--index_file', dest="index_file", required=True, help="index file with baits information")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()

