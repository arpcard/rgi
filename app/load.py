import shutil
import argparse
from app.settings import *

# this script is used to load new card.json file to system wide package or local

def validate_file(filename):
	try:
		with open(filename) as f:
			out = json.load(f)
			return out
	except ValueError as e:
		logger.error("Invalid json: {}".format(e))
		return None # or: raise

def main(args):
	# print args
	logger.info(json.dumps(args.__dict__, indent=2))
	if args.debug:
		logger.setLevel(10)

	if args.card_json is not None:
		# validate json
		if validate_file(args.card_json) == False:
			logger.error("failed to read json file: {}".format(args.card_json))
			exit()
		load_file(args.local_database, args.card_json, "card.json")

	if args.card_annotation is not None and (args.wildcard_index is None or args.wildcard_annotation is None):
		load_reference_card_only(args.local_database, args.card_annotation, "card_reference.fasta")

	if args.wildcard_index is not None and args.wildcard_annotation is not None and args.card_annotation is not None:
		# load index
		load_file(args.local_database, args.wildcard_index, "index-for-model-sequences.txt")
		# load annotation files (card and wildcard)
		load_reference_card_and_wildcard(args.local_database, args.card_annotation , args.wildcard_annotation,"card_wildcard_reference.fasta")

	if args.kmer_database is not None:
		if args.kmer_size is not None:
			load_file(args.local_database, args.kmer_database, "{}mer_database.json".format(str(args.kmer_size)))
		else:
			logger.error("Need to specify kmer size when loading kmer files.")

	if args.amr_kmers is not None:
		if args.kmer_size is not None:
			load_file(args.local_database, args.amr_kmers, "amr_{}mer.txt".format(str(args.kmer_size)))
		else:
			logger.error("Need to specify kmer size when loading kmer files.")

def load_reference_card_only(local_db, fasta_file, filename):
	load_file(local_db, fasta_file, "card_reference.fasta")
	logger.info("loaded card only for 'rgi bwt'.")

def load_reference_card_and_wildcard(local_db, card_fasta_file, wildcard_fasta_file, filename):

	db = get_location(local_db)

	filenames = []
	filenames.append(card_fasta_file)
	filenames.append(wildcard_fasta_file)

	# combine the two fastas
	import fileinput
	with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
		for line in fin:
			fout.write(line)
	logger.info("loaded card and wildcard annotations for 'rgi bwt'.")

def get_location(local_db):
	db = ""

	# path to save card.json file and database files
	if local_db == True:
		db = LOCAL_DATABASE
		# create directory if it doesn't exist
		if not os.path.exists(LOCAL_DATABASE):
			os.makedirs(LOCAL_DATABASE)
	else:
		db = data_path

	return db

def load_file(local_db, filepath, filename, validate_json=False):

	db = get_location(local_db)

	try:
		# copy new file
		shutil.copyfile(filepath, os.path.join(db, filename))
		logger.info("file {} loaded ok".format(filename))
	except Exception as e:
		logger.warning("failed to copy json file: {}".format(e))

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('-i', '--card_json', required=False, help='must be a card database json file')

	parser.add_argument('--card_annotation', required=False, help="annotated reference FASTA")

	parser.add_argument('--wildcard_annotation', required=False, help="annotated reference FASTA")
	parser.add_argument('--wildcard_index', required=False, help="wildcard index file (index-for-model-sequences.txt)")

	parser.add_argument('--kmer_database', required=False, help="json of kmer database")
	parser.add_argument('--amr_kmers', required=False, help="txt file of all amr kmers")
	parser.add_argument('--kmer_size', required=False, help="kmer size if loading kmer files")

	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()
