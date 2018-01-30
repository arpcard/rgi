from app.settings import *
from app.RGI import RGI
import argparse
from app.ConvertRGIJsonToTSV import ConvertJsonToTSV
from app.Galaxy import Galaxy
import app.Parser
import app.load
import app.clean

class MainBase(object):
    def __init__(self, api=False):
        # """
        USAGE='''%(prog)s <command> [<args>] 
            commands are:
               main    Runs rgi application
               tab     Creates a Tab-delimited from rgi results
               parser  Creates categorical .json files RGI wheel visualization. An input .json file containing the RGI results must be input.
               load    Loads CARD database json file
               clean   Removes BLAST databases and temporary files
               galaxy  Galaxy project wrapper'''

        parser = argparse.ArgumentParser(prog="rgi", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        parser.add_argument('command', help='Subcommand to run')

        if api == False:
            args=parser.parse_args(sys.argv[1:2])
            # """
            if not hasattr(self, args.command):
                logger.info("Unrecognized command: {}".format(args.command))
                # parser.print_help()
                exit("Error: Unrecognized command: {}".format(args.command))
            getattr(self, args.command)()

    def main(self):
        parser = self.main_args()
        args = parser.parse_args(sys.argv[2:])
        self.main_run(args)

    def main_args(self):
        parser = argparse.ArgumentParser(prog="rgi main",description="{} - {} - Main".format(APP_NAME,SOFTWARE_VERSION))
        parser.add_argument('-i','--input_sequence', dest="input_sequence", default=None, required=True, \
            help='input file must be in either FASTA (contig and protein), FASTQ(read) or gzip format! e.g myFile.fasta, myFasta.fasta.gz')
        parser.add_argument('-o','--output_file', dest="output_file", default=None,required=True, help="output JSON file (default=None)")
        parser.add_argument('-t','--input_type', dest="input_type", default="CONTIG", required=False, help='must be one of contig, protein, read (default: contig)')
        parser.add_argument('-a','--alignment_tool', dest="aligner", default="BLAST", help = "specify alignment tool. Options are BLAST or DIAMOND  (default = BLAST)")
        parser.add_argument('-n','--num_threads', dest="threads", default="32", help="number of threads (CPUs) to use in the BLAST search (default=32)")
        parser.add_argument('--include_loose', dest="loose", action='store_true', help="include loose hits in addition to strict and perfect hits")
        parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files")
        parser.add_argument('-d','--data', dest="data", default="NA", help = "specify a data-type, i.e. wgs, chromosome, plasmid, etc. (default = NA)")
        parser.add_argument('-v','--software_version', action='version', version="{}".format(SOFTWARE_VERSION), help = "prints software number")
        parser.add_argument('-dv','--data_version',action='version',version="{}".format(DATA_VERSION), help = "prints data version number")
        return parser

    def main_run(self, args):
        logger.info('Running rgi main')
        rgi_obj = RGI(**vars(args))
        rgi_obj.run()

    def tab(self):
        parser = self.tab_args()
        args = parser.parse_args(sys.argv[2:])
        self.tab_run(args)

    def tab_args(self):
        parser = argparse.ArgumentParser(prog="rgi tab", description="{} - {} - Tab-delimited".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-i', '--afile', help='must be a rgi json result file', required=True)
        return parser

    def tab_run(self, args):
        logger.info('Running rgi tab')
        obj = ConvertJsonToTSV(args.afile)
        obj.run()

    def parser(self):
        parser = self.parser_args()
        args = parser.parse_args(sys.argv[2:])
        self.parser_run(args)

    def parser_args(self):
        parser = app.Parser.create_parser()
        return parser

    def parser_run(self, args):
        app.Parser.api_main(args)

    def load(self):
        parser = self.load_args()
        args = parser.parse_args(sys.argv[2:])
        self.load_run(args)

    def load_args(self):
        parser = argparse.ArgumentParser(prog="rgi load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-i','--afile', help='must be a card database json file', required=True)
        return parser

    def load_run(self, args):
        logger.info('Running rgi load')
        app.load.main(args)

    def clean(self):
        parser = self.clean_args()
        args = parser.parse_args(sys.argv[2:])
        self.clean_run(args)

    def clean_args(self):
        parser = argparse.ArgumentParser(prog="rgi clean", description="{} - {} - Clean".format(APP_NAME, SOFTWARE_VERSION))
        return parser

    def clean_run(self, args):
        logger.info('Running rgi clean')
        app.clean.main()

    def galaxy(self):
        parser = self.galaxy_args()
        args = parser.parse_args(sys.argv[2:])
        self.galaxy_run(args)

    def galaxy_args(self):
        parser = argparse.ArgumentParser(prog="rgi galaxy", description="{} - {} - Galaxy project wrapper".\
            format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('--galaxy_database', help='specify path to CARD blast databases. Note this is \
            ONLY used for galaxyproject wrapper (default: None)', required=True)
        return parser

    def galaxy_run(self, args):
        logger.info('Running rgi galaxy')
        obj = Galaxy(args.galaxy_database)
        obj.load_db_galaxy()


if __name__ == '__main__':
    MainBase()

