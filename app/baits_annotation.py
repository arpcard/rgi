import os, sys, json, csv, argparse
from app.settings import *

def main(args):
    working_directory = os.getcwd()
    print(args)

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

