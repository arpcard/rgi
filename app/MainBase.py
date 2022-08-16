from app.settings import *
from app.RGI import RGI
import argparse
from app.ConvertRGIJsonToTSV import ConvertJsonToTSV
from app.Galaxy import Galaxy
import app.Parser
import app.load
import app.auto_load
import app.clean
import app.build_kmer_sets
import app.card_annotation
import app.wildcard_annotation
import app.baits_annotation
import app.remove_duplicates
from app.kmer_query import CARDkmers
from app.BWT import BWT
from app.Heatmap import Heatmap
from app.Baits import Baits
from argparse import RawTextHelpFormatter

class MainBase(object):
    def __init__(self, api=False):
        # """
        self.cpu_count = os.cpu_count()
        USAGE='''%(prog)s <command> [<args>]
            commands are:
               ---------------------------------------------------------------------------------------
               Database
               ---------------------------------------------------------------------------------------
               auto_load Automatically loads CARD database, annotations and k-mer database
               load      Loads CARD database, annotations and k-mer database
               clean     Removes BLAST databases and temporary files
               database  Information on installed card database
               galaxy    Galaxy project wrapper

               ---------------------------------------------------------------------------------------
               Genomic
               ---------------------------------------------------------------------------------------

               main     Runs rgi application
               tab      Creates a Tab-delimited from rgi results
               parser   Creates categorical JSON files RGI wheel visualization
               heatmap  Heatmap for multiple analysis

               ---------------------------------------------------------------------------------------
               Metagenomic
               ---------------------------------------------------------------------------------------
               bwt                   Align reads to CARD and in silico predicted allelic variants (beta)

               ---------------------------------------------------------------------------------------
               Baits validation
               ---------------------------------------------------------------------------------------
               tm                    Baits Melting Temperature

               ---------------------------------------------------------------------------------------
               Annotations
               ---------------------------------------------------------------------------------------
               card_annotation       Create fasta files with annotations from card.json
               wildcard_annotation   Create fasta files with annotations from variants
               baits_annotation      Create fasta files with annotations from baits (experimental)
               remove_duplicates     Removes duplicate sequences (experimental)


               ---------------------------------------------------------------------------------------
               Pathogen of origin
               ---------------------------------------------------------------------------------------

               kmer_build            Build AMR specific k-mers database used for pathogen of origin (beta)
               kmer_query            Query sequences against AMR k-mers database to predict pathogen of origin (beta)

               '''

        parser = argparse.ArgumentParser(prog="rgi", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        parser.add_argument('command', choices=['main', 'tab', 'parser', 'load', 'auto_load',
                                                'clean', 'galaxy', 'database', 'bwt', 'tm', 'card_annotation', 'wildcard_annotation', 'baits_annotation', 'remove_duplicates', 'heatmap', 'kmer_build', 'kmer_query'],
                                                help='Subcommand to run')

        if api == False:
            args=parser.parse_args(sys.argv[1:2])
            if not hasattr(self, args.command):
                logger.info("Unrecognized command: {}".format(args.command))
                exit("Error: Unrecognized command: {}".format(args.command))
            getattr(self, args.command)()

    def main(self):
        parser = self.main_args()
        args = parser.parse_args(sys.argv[2:])
        self.main_run(args)

    def main_args(self):
        parser = argparse.ArgumentParser(prog="rgi main",description="{} - {} - Main".format(APP_NAME,SOFTWARE_VERSION))
        parser.add_argument('-i','--input_sequence', dest="input_sequence", required=True, \
            help='input file must be in either FASTA (contig and protein) or gzip format! e.g myFile.fasta, myFasta.fasta.gz')
        parser.add_argument('-o','--output_file', dest="output_file", required=True, help="output folder and base filename")
        parser.add_argument('-t','--input_type', dest="input_type",
                type=str.lower,
                default="contig", choices=['contig','protein'],
                required=False,
                help='specify data input type (default = contig)')
        parser.add_argument('-a','--alignment_tool', dest="aligner",
                type=str.upper,
                choices = ['DIAMOND', 'BLAST'],
                default="BLAST",
                help = "specify alignment tool (default = BLAST)")
        parser.add_argument('-n','--num_threads', dest="threads", type=int,
                default=self.cpu_count, help="number of threads (CPUs) to use in the BLAST search (default={})".format(self.cpu_count))
        parser.add_argument('--include_loose', dest="loose", action='store_true', help="include loose hits in addition to strict and perfect hits (default: False)")
        parser.add_argument('--include_nudge', dest="include_nudge", action='store_true', help="include hits nudged from loose to strict hits (default: False)")
        parser.add_argument('--local', dest="local_database", action='store_true', help="use local database (default: uses database in executable directory)")
        parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files (default: False)")
        parser.add_argument('--keep', dest="keep", action="store_true", help="keeps Prodigal CDS when used with --clean (default: False)")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode (default: False)")
        parser.add_argument('--low_quality', dest="low_quality", action="store_true", help="use for short contigs to predict partial genes (default: False)")
        parser.add_argument('-d','--data', dest="data", default="NA",
                choices=['wgs', 'plasmid', 'chromosome', 'NA'],
                help = "specify a data-type (default = NA)")
        parser.add_argument('-v','--version', action='version', version="{}".format(SOFTWARE_VERSION), help = "prints software version number")
        parser.add_argument('--split_prodigal_jobs', dest="split_prodigal_jobs", action="store_true", help="run multiple prodigal jobs simultaneously for contigs in a fasta file (default: False)")
        return parser

    def main_run(self, args):
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
        parser = app.load.create_parser()
        return parser

    def load_run(self, args):
        app.load.main(args)

    def auto_load(self):
        parser = self.auto_load_args()
        args = parser.parse_args(sys.argv[2:])
        self.auto_load_run(args)

    def auto_load_args(self):
        parser = app.auto_load.create_parser()
        return parser

    def auto_load_run(self, args):
        app.auto_load.main(args)

    def kmer_build(self):
        parser = self.kmer_build_args()
        args = parser.parse_args(sys.argv[2:])
        self.kmer_build_run(args)

    def kmer_build_args(self):
        parser = app.build_kmer_sets.create_parser()
        return parser

    def kmer_build_run(self, args):
        app.build_kmer_sets.main(args)

    def kmer_query(self):
        parser = self.kmer_query_args()
        args = parser.parse_args(sys.argv[2:])
        self.kmer_query_run(args)

    def kmer_query_args(self):
        parser = argparse.ArgumentParser(prog="rgi kmer_query",
            description='{} - {} - Kmer Query \n\nTests sequenes using CARD*kmers'.format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
        parser.add_argument('-i', '--input', dest="input", required=True,
            help="Input file (bam file from RGI*BWT, json file of RGI results, fasta file of sequences)")
        parser.add_argument('--bwt', action="store_true",
            help="Specify if the input file for analysis is a bam file generated from RGI*BWT")
        parser.add_argument('--rgi', action="store_true",
            help="Specify if the input file is a RGI results json file")
        parser.add_argument('--fasta', action="store_true",
            help="Specify if the input file is a fasta file of sequences")
        parser.add_argument('-k', '--kmer_size', dest="k", required=True,
            help="length of k")
        parser.add_argument('-m', '--minimum', dest="min", default=10,
            help="Minimum number of kmers in the called category for the classification to be made (default=10).")
        parser.add_argument('-n','--threads', dest="threads", type=int,
            default=1, help="number of threads (CPUs) to use (default={})".format(1))
        parser.add_argument('-o', '--output', dest="output", required=True,
            help="Output file name.")
        parser.add_argument('--local', dest="local_database", action='store_true',
            help="use local database (default: uses database in executable directory)")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
        return parser

    def kmer_query_run(self, args):
        obj = CARDkmers(args.input, args.bwt, args.rgi, args.fasta, args.k, args.min, args.threads, args.output, args.local_database, args.debug)
        obj.run()

    def card_annotation(self):
        parser = self.card_annotation_args()
        args = parser.parse_args(sys.argv[2:])
        self.card_annotation_run(args)

    def card_annotation_args(self):
        parser = app.card_annotation.create_parser()
        return parser

    def card_annotation_run(self, args):
        app.card_annotation.main(args)

    def wildcard_annotation(self):
        parser = self.wildcard_annotation_args()
        args = parser.parse_args(sys.argv[2:])
        self.wildcard_annotation_run(args)

    def wildcard_annotation_args(self):
        parser = app.wildcard_annotation.create_parser()
        return parser

    def wildcard_annotation_run(self, args):
        app.wildcard_annotation.main(args)

    def baits_annotation(self):
        parser = self.baits_annotation_args()
        args = parser.parse_args(sys.argv[2:])
        self.baits_annotation_run(args)

    def baits_annotation_args(self):
        parser = app.baits_annotation.create_parser()
        return parser

    def baits_annotation_run(self, args):
        app.baits_annotation.main(args)

    def remove_duplicates(self):
        parser = self.remove_duplicates_args()
        args = parser.parse_args(sys.argv[2:])
        self.remove_duplicates_run(args)

    def remove_duplicates_args(self):
        parser = app.remove_duplicates.create_parser()
        return parser

    def remove_duplicates_run(self, args):
        app.remove_duplicates.main(args)

    def bwt(self):
        parser = self.bwt_args()
        args = parser.parse_args(sys.argv[2:])
        self.bwt_run(args)

    def bwt_args(self):
		# description="{} - {} - Main".format(APP_NAME,SOFTWARE_VERSION))
        parser = argparse.ArgumentParser(prog="rgi bwt",description="{} - {} - BWT \n\nAligns metagenomic reads to CARD and wildCARD reference using kma, bowtie2 or bwa and provide reports.".format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
        parser.add_argument('-1', '--read_one', required=True, help="raw read one (qc and trimmed)")
        parser.add_argument('-2', '--read_two', help="raw read two (qc and trimmed)")
        parser.add_argument('-a', '--aligner', default="kma", choices=['kma','bowtie2','bwa'], help="select read aligner (default=kma)")
        parser.add_argument('-n','--threads', dest="threads", type=int,default=self.cpu_count, help="number of threads (CPUs) to use (default={})".format(self.cpu_count))
        parser.add_argument('-o','--output_file', dest="output_file", required=True, help="name of output filename(s)")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode (default=False)")
        parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files (default=False)")
        parser.add_argument('--local', dest="local_database", action='store_true', help="use local database (default: uses database in executable directory)")
        parser.add_argument('--include_wildcard', dest="include_wildcard", action="store_true", help="include wildcard (default=False)")
        parser.add_argument('--include_other_models', dest="include_other_models", action="store_true", \
            help="include protein variant, rRNA variant, knockout, and protein overexpression models (default=False)")
        parser.add_argument('--include_baits', dest="include_baits", action="store_true", help="include baits (default=False)")
        parser.add_argument('--mapq', dest="mapq", help="filter reads based on MAPQ score (default=False)")
        parser.add_argument('--mapped', dest="mapped", help="filter reads based on mapped reads (default=False)")
        parser.add_argument('--coverage', dest="coverage", help="filter reads based on coverage of reference sequence")

        return parser

    def bwt_run(self, args):
        obj = BWT(
            args.aligner,
            args.include_wildcard,
            args.include_baits,
            args.read_one,
            args.read_two,
            args.threads,
            args.output_file,
            args.debug,
            args.clean,
            args.local_database,
            args.mapq,
            args.mapped,
            args.coverage,
            args.include_other_models
        )
        obj.run()

    def tm(self):
        parser = self.tm_args()
        args = parser.parse_args(sys.argv[2:])
        self.tm_run(args)

    def tm_args(self):
        parser = argparse.ArgumentParser(prog="rgi tm",description='{} - {} - TM'.format(APP_NAME,SOFTWARE_VERSION))
        parser.add_argument('-i', '--input_file', dest="input_file", help="input_file")
        parser.add_argument('-o', '--output_file', dest="output_file", help="output_file")
        parser.add_argument('-t', '--filter_temperature', dest="filter_temperature", default=65,
            help="desired melting temperature (default=65).")
        parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
        return parser

    def tm_run(self, args):
        obj = Baits(
            args.input_file,
            args.output_file,
            args.filter_temperature,
            args.clean,
            args.debug
        )
        obj.run()

    def heatmap(self):
        parser = self.heatmap_args()
        args = parser.parse_args(sys.argv[2:])
        self.heatmap_run(args)

    def heatmap_args(self):
        parser = argparse.ArgumentParser(prog="rgi heatmap",description='{} - {} - Heatmap \n\nCreates a heatmap when given multiple RGI results.'.format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
        parser.add_argument('-i', '--input', dest="input", required=True, help="Directory containing the RGI .json files (REQUIRED)")
        parser.add_argument('-cat', '--category', dest="classification", choices=("drug_class", "resistance_mechanism", "gene_family"),
            help="The option to organize resistance genes based on a category.")
        parser.add_argument('-f', '--frequency', dest="frequency", action="store_true", help="Represent samples based on resistance profile.")
        parser.add_argument('-o', '--output', dest="output", default="RGI_heatmap", help="Name for the output EPS and PNG files.\nThe number of files run will automatically \nbe appended to the end of the file name.(default={})".format('RGI_heatmap'))
        parser.add_argument('-clus', '--cluster', dest="cluster", choices=("samples", "genes", "both"),
            help="Option to use SciPy's hiearchical clustering algorithm to cluster rows (AMR genes) or columns (samples).")
        parser.add_argument('-d', '--display', dest="display", choices=("plain", "fill", "text"), default="plain",
            help="Specify display options for categories (deafult=plain).")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")

        return parser

    def heatmap_run(self, args):
        obj = Heatmap(args.input, args.classification, args.frequency, args.output, args.cluster, args.display, args.debug)
        obj.run()

    def clean(self):
        parser = self.clean_args()
        args = parser.parse_args(sys.argv[2:])
        self.clean_run(args)

    def clean_args(self):
        parser = app.clean.create_parser()
        return parser

    def clean_run(self, args):
        app.clean.main(args)

    def galaxy(self):
        parser = self.galaxy_args()
        args = parser.parse_args(sys.argv[2:])
        self.galaxy_run(args)

    def galaxy_args(self):
        parser = argparse.ArgumentParser(prog="rgi galaxy", description="{} - {} - Galaxy project wrapper".\
            format(APP_NAME, SOFTWARE_VERSION), epilog=GALAXY_PROJECT_WRAPPER)
        parser.add_argument('--galaxy_database', help='path to CARD data and blast databases', required=True)
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
        return parser

    def galaxy_run(self, args):
        obj = Galaxy(args.galaxy_database, args.debug)
        obj.load_db_galaxy()

    def database(self):
        parser = self.database_args()
        args = parser.parse_args(sys.argv[2:])
        print(self.database_run(args))

    def database_args(self):
        parser = argparse.ArgumentParser(prog="rgi database", description="{} - {} - Database".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-v','--version',action='store_true', required=True, help = "prints data version number")
        parser.add_argument('--local', dest="local_database", action='store_true', help="use local database (default: uses database in executable directory)")
        parser.add_argument('--all', action='store_true', help="data version number used for `rgi bwt` and `rgi main` (default: rgi main)")
        return parser

    def database_run(self, args):
        data_version = ""
        # path to loaded database files
        if args.local_database:
            db = LOCAL_DATABASE
            # error if it doesn't exist
            if not os.path.exists(LOCAL_DATABASE):
                print("Error: missing local directory: {}".format(LOCAL_DATABASE))
                print("Please run `rgi load --local -i <path to card.json>` to create local database.")
                print("See `rgi load --help` to upload the card.json to rgi application.\n".format(os.path.abspath(LOCAL_DATABASE)))
                exit()
        else:
            db = data_path

        indecies_directory = os.path.join(db,"loaded_databases.json")

        if os.path.isfile(indecies_directory) == True:
            with open(indecies_directory) as json_file:
                json_data = json.load(json_file)
                if args.all == True:
                    kmers_str = ",".join(json_data["card_kmers"]["kmer_sizes"])
                    if kmers_str == "":
                        kmers_str = "N/A"
                    data_version = ("card_canonical: {} | card_variants: {} | kmer_sizes: {}".format(
                        json_data["card_canonical"]["data_version"],
                        json_data["card_variants"]["data_version"],
                        kmers_str
                        )
                    )
                    if "model_type_used" in json_data["card_canonical"].keys() or "model_type_used" in json_data["card_variants"].keys() :
                        card_canonical_model_type_used = "N/A"
                        card_variants_model_type_used = "N/A"

                        if len(json_data["card_canonical"]["model_type_used"]) > 0:
                            card_canonical_model_type_used = ";".join(json_data["card_canonical"]["model_type_used"])
                        if len(json_data["card_variants"]["model_type_used"]) > 0:
                            card_variants_model_type_used = ";".join(json_data["card_variants"]["model_type_used"])

                        data_version = ("card_canonical: {} | card_canonical_model_type_used: {} | card_variants: {} | card_variants_model_type_used: {} | kmer_sizes: {}".format(
                            json_data["card_canonical"]["data_version"],
                            card_canonical_model_type_used,
                            json_data["card_variants"]["data_version"],
                            card_variants_model_type_used,
                            kmers_str
                            )
                        )
                else:
                    data_version = json_data["card_canonical"]["data_version"]
        else:
            print('\nError: no databases found in data path: {}. \nSee `rgi load --help`\n'.format(os.path.abspath(db)))
            exit()

        if data_version == "":
            print('\nError: no databases found in data path: {}. \nSee `rgi load --help`\n'.format(os.path.abspath(db)))

        return data_version

if __name__ == '__main__':
    MainBase()
