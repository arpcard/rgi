import shutil
import argparse
from app.settings import *
import tempfile
from bs4 import BeautifulSoup as bs
import requests, json
import re

def main(args):
	debug = ""
	if args.debug:
		logger.setLevel(10)
		debug = "--debug"

	local_database = ""
	if args.local_database:
		local_database = "--local"

	card_cannonical_version, card_variants_version = get_versions()

	logger.info(json.dumps(args.__dict__, indent=2))
	logger.info("card cannonical version: {}".format(card_cannonical_version))
	logger.info("card variants version: {}".format(card_variants_version))

	# Create the directory
	logger.info("Creating temp directory with prefix {}".format(os.path.join(os.getcwd(), "rgi_autoload_")))
	directory = tempfile.mkdtemp(prefix=os.path.join(os.getcwd(), "rgi_autoload_"))
	logger.info("Directory '%s' created" % directory)
	# print("Directory '%s' created" % directory)
	logger.info("=================================== DOWNLOAD CARD CANONICAL DATA ===================================")
	# get latest card database
	data=os.path.join(directory,"data")
	card_data=os.path.join(directory,"card_data")
	os.system("wget -O {data} --no-check-certificate https://card.mcmaster.ca/download/0/broadstreet-v{card_cannonical_version}.tar.bz2".format(
		data=data,
		card_cannonical_version=card_cannonical_version
		)
	)
	os.system("mkdir -p {card_data}".format(card_data=card_data))
	os.system("tar xf {data} -C {card_data}".format(data=data,card_data=card_data))

	logger.info("=================================== DOWNLOAD CARD VARIANTS DATA ===================================")
	variants=os.path.join(directory,"variants")
	card_variants=os.path.join(directory,"card_variants")
	os.system("wget -O {variants} --no-check-certificate https://card.mcmaster.ca/download/6/prevalence-v{card_variants_version}.tar.bz2".format(
		variants=variants,
		card_variants_version=card_variants_version
		)
	)
	os.system("mkdir -p {card_variants}".format(card_variants=card_variants))
	os.system("tar xf {variants} -C {card_variants}".format(variants=variants,card_variants=card_variants))
	os.system("gunzip {card_variants}/*.gz".format(card_variants=card_variants))

	logger.info("=================================== CARD CANONICAL ANNOTATIONS ===================================")
	os.system("rgi card_annotation --input {card_data}/card.json".format(card_data=card_data))

	logger.info("=================================== CARD VARIANTS ANNOTATIONS ===================================")
	os.system("rgi wildcard_annotation --input_directory {card_variants} --version {card_variants_version} --card_json {card_data}/card.json".format(
		card_variants=card_variants,
		card_variants_version=card_variants_version,
		card_data=card_data
		)
	)

	logger.info("=================================== CLEAN OLD DATABASES ===================================")
	os.system("rgi clean {debug} {local_database}".format(local_database=local_database,debug=debug))

	logger.info("=================================== LOAD DATABASES ===================================")
	os.system("rgi load \
	--card_json {card_data}/card.json \
	--card_annotation card_database_v{card_cannonical_version}.fasta \
	--card_annotation_all_models card_database_v{card_cannonical_version}_all.fasta \
	--wildcard_index {card_variants}/index-for-model-sequences.txt \
	--wildcard_version {card_variants_version} \
	--wildcard_annotation wildcard_database_v{card_variants_version}.fasta \
	--wildcard_annotation_all_models wildcard_database_v{card_variants_version}_all.fasta \
	--kmer_database {card_variants}/61_kmer_db.json \
	--amr_kmers {card_variants}/all_amr_61mers.txt \
	--kmer_size 61 \
	{local_database} {debug}".format(
		card_data=card_data,
		card_cannonical_version=card_cannonical_version,
		card_variants=card_variants,
		card_variants_version=card_variants_version,
		local_database=local_database,
		debug=debug
		)
	)

	logger.info("=================================== CHECK LOADED DATABASES ===================================")
	os.system("rgi database -v --all {local_database}".format(local_database=local_database))

	if args.clean:
		logger.info("=================================== CLEAN UP ===================================")
		os.system("rm {data}".format(data=data))
		os.system("rm {variants}".format(variants=variants))
		os.system("rm {card_data}/* ".format(card_data=card_data))
		os.system("rm {card_variants}/* ".format(card_variants=card_variants))
		os.system("rm -r {card_data} ".format(card_data=card_data))
		os.system("rm -r {card_variants}".format(card_variants=card_variants))
		os.system("rm -r {}".format(directory))
		os.system("rm card_database_v{card_cannonical_version}.fasta".format(card_cannonical_version=card_cannonical_version))
		os.system("rm card_database_v{card_cannonical_version}_all.fasta".format(card_cannonical_version=card_cannonical_version))
		os.system("rm wildcard_database_v{card_variants_version}.fasta".format(card_variants_version=card_variants_version))
		os.system("rm wildcard_database_v{card_variants_version}_all.fasta".format(card_variants_version=card_variants_version))

	logger.info("=================================== DONE ===================================")

def get_versions():
	r = requests.get('https://card.mcmaster.ca/download')
	soup = bs(r.content, 'html.parser')
	data = [item['href'] if item.get('href') is not None else item['src'] for item in soup.select('[href^="/download/0"]') ]
	data_version = valid_version(re.search(r'\s*([\d.].([\d.]).([\d.]))', data[0]).group(1))
	prev = [item['href'] if item.get('href') is not None else item['src'] for item in soup.select('[href^="/download/6"]') ]
	prev_version = valid_version(re.search(r'\s*([\d.].([\d.]).([\d.]))', prev[0]).group(1))
	return (data_version, prev_version)

def valid_version(s):
	msg = "Not a version: '{0}'.".format(s)
	pattern = re.compile('^\s*([\d.]).([\d.]).([\d.])$')
	if pattern.match(s) is not None:
		return s
	else:
		raise argparse.ArgumentTypeError(msg)

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi auto_load", description="{} - {} - Automatic Load".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()
