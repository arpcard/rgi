import shutil
import argparse
from app.settings import *

'''
loaded_databases example:

{
    "card_json": {
      "data_version": "3.0.1",
	  "model_type_used": [
	      "homolog",
	      "variant",
	      "rRNA",
	      "overexpression",
	      "knockout"
	    ]
    },
    "card_canonical": {
      "data_version": "3.0.1",
	  "model_type_used": [
	      "homolog",
	      "variant",
	      "rRNA",
	      "overexpression",
	      "knockout"
	    ]
    },
    "card_variants": {
      "data_version": "4.0.0",
	  "model_type_used": [
	      "homolog",
	      "variant",
	      "rRNA",
	      "overexpression",
	      "knockout"
	    ]
    },
    "card_kmers": {
      "kmer_sizes": ["31","61"]
    }
  }
'''
loaded_databases = {}

# this script is used to load new card.json file to system wide package or local
def get_card_json_version(card_json_path):
	data_version = ""
	if os.path.isfile(card_json_path) == True:
		with open(card_json_path) as json_file:
			json_data = json.load(json_file)
			for item in json_data.keys():
				if item == "_version":
					data_version = json_data[item]
	return data_version

def validate_file(filename):
	try:
		with open(filename) as f:
			out = json.load(f)
			return out
	except ValueError as e:
		logger.error("Invalid json: {}".format(e))
		return None # or: raise

def get_card_annotation(args_card_annotation, args_local_database, filename="card_reference.fasta"):
	card_annotation = args_card_annotation
	# check if card annotation is already loaded
	if args_card_annotation is None:
		db = get_location(args_local_database)
		if os.path.isfile(os.path.join(db, filename)):
			logger.info("card annotation file {} exists".format(filename))
			card_annotation = os.path.join(db, filename)
	return card_annotation

def main(args):
	# versions used
	# print args
	if args.debug:
		logger.setLevel(10)

	version_file = os.path.join(get_location(args.local_database), "loaded_databases.json")
	if os.path.isfile(version_file):
		# file exists load
		logger.info("file {} exists load".format(version_file))
		with open(version_file, 'r') as fout:
			loaded_databases = json.load(fout)
	else:
		# initialize bwt dict for versions used
		loaded_databases = {
				"card_json": {
					"data_version": "N/A",
					"model_type_used": []
				},
				"card_canonical": {
					"data_version": "N/A",
					"model_type_used": []
				},
				"card_variants": {
					"data_version": "N/A",
					"model_type_used": []
				},
				"card_kmers": {
					"kmer_sizes": []
				}
			}

	logger.info(json.dumps(args.__dict__, indent=2))

	# Load card.json used for `rgi main` which uses all models by default
	if args.card_json is not None:
		# validate json
		if validate_file(args.card_json) == False:
			logger.error("failed to read json file: {}".format(args.card_json))
			exit()
		load_file(args.local_database, args.card_json, "card.json")

		# ['card_json', 'card_canonical', 'card_variants', 'card_kmers']
		loaded_databases["card_json"]["data_version"] = get_card_json_version(os.path.join(get_location(args.local_database), "card.json"))
		loaded_databases["card_json"]["model_type_used"] = ["homolog","variant","rRNA","overexpression","knockout"]

	card_annotation = args.card_annotation
	card_annotation_all = args.card_annotation_all_models

	# if args.card_annotation is not None and (args.wildcard_index is None or args.wildcard_annotation is None):
	# CARD[homolog]
	if args.card_annotation is not None:
		load_reference_card_only(args.local_database, args.card_annotation, "card_reference.fasta")
		loaded_databases["card_canonical"]["model_type_used"] = ["homolog"]
		loaded_databases["card_canonical"]["data_version"] = loaded_databases["card_json"]["data_version"]
	else:
		# check if card annotation is already loaded
		card_annotation = get_card_annotation(args.card_annotation, args.local_database)

	# CARD[All]
	if args.card_annotation_all_models is not None:
		load_reference_card_only(args.local_database, args.card_annotation_all_models, "card_reference_all.fasta")
		loaded_databases["card_canonical"]["model_type_used"] = ["homolog","variant","rRNA","overexpression","knockout"]
		loaded_databases["card_canonical"]["data_version"] = loaded_databases["card_json"]["data_version"]
	else:
		card_annotation_all = get_card_annotation(args.card_annotation_all_models, args.local_database, "card_reference_all.fasta")

	# load wildcard index file
	if args.wildcard_index is not None:
		load_file(args.local_database, args.wildcard_index, "index-for-model-sequences.txt")

	# CARD[homolog] + WILD[homolog]
	if args.wildcard_annotation is not None and args.card_annotation is not None:
		# load annotation files (card and wildcard)
		load_reference_card_and_wildcard(args.local_database, card_annotation , args.wildcard_annotation,"card_wildcard_reference.fasta")
		loaded_databases["card_variants"]["data_version"] = args.wildcard_version
		loaded_databases["card_variants"]["model_type_used"] = ["homolog"]

	# CARD[All] + WILD[ALL]
	if args.wildcard_annotation_all_models is not None and args.card_annotation_all_models is not None:
		load_reference_card_and_wildcard(args.local_database, card_annotation_all , args.wildcard_annotation_all_models,"card_wildcard_reference_all.fasta")
		loaded_databases["card_variants"]["data_version"] = args.wildcard_version
		loaded_databases["card_variants"]["model_type_used"] = ["homolog","variant","rRNA","overexpression","knockout"]

	if args.kmer_database is not None:
		if args.kmer_size is not None:
			load_file(args.local_database, args.kmer_database, "{}mer_database.json".format(str(args.kmer_size)))
			if args.kmer_size not in loaded_databases["card_kmers"]["kmer_sizes"]:
				loaded_databases["card_kmers"]["kmer_sizes"].append(args.kmer_size)
		else:
			logger.error("Need to specify kmer size when loading kmer files.")

	if args.amr_kmers is not None:
		if args.kmer_size is not None:
			load_file(args.local_database, args.amr_kmers, "amr_{}mer.txt".format(str(args.kmer_size)))
			if args.kmer_size not in loaded_databases["card_kmers"]["kmer_sizes"]:
				loaded_databases["card_kmers"]["kmer_sizes"].append(args.kmer_size)
		else:
			logger.error("Need to specify kmer size when loading kmer files.")

	if args.baits_index is not None and args.baits_annotation is not None:#and args.card_annotation is not None:
		logger.info("adding index and fasta for baits")
		# load index
		load_file(args.local_database, args.baits_index, "baits-probes-with-sequence-info.txt")
		# load annotation files (baits)
		load_reference_card_and_baits(args.local_database, card_annotation, args.baits_annotation, "card_baits_reference.fasta")

	if args.wildcard_index is not None and args.wildcard_annotation is not None and args.baits_index is not None and args.baits_annotation is not None:
		# load annotations files for CARD, VARIANTS and BAITS
		load_reference_card_and_wilcard_and_baits(args.local_database, card_annotation, args.baits_annotation, args.wildcard_annotation,  "card_wildcard_baits_reference.fasta")
		loaded_databases["card_variants"]["data_version"] = args.wildcard_version

	# write versions used
	with open(version_file, 'w') as fout:
		json.dump(loaded_databases, fout)
	# print out loaded databases
	logger.info(json.dumps(loaded_databases, indent=2))

def load_reference_card_only(local_db, fasta_file, filename):
	# load_file(local_db, fasta_file, "card_reference.fasta")
	load_file(local_db, fasta_file, filename)
	logger.info("loaded card annotation file: {} for 'rgi bwt'.".format(filename))

def load_reference_baits_only(local_db, fasta_file):
	load_file(local_db, fasta_file, "baits_reference.fasta")
	logger.info("loaded baits only for 'rgi bwt'.")

def load_reference_card_and_wilcard_and_baits(local_db, card_fasta_file, baits_fasta_file, wildcard_fasta_file, filename):
	db = get_location(local_db)
	filenames = []
	filenames.append(card_fasta_file)
	filenames.append(wildcard_fasta_file)
	filenames.append(baits_fasta_file)

	# combine the three fastas
	import fileinput
	with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
		for line in fin:
			fout.write(line)
	logger.info("loaded card, wildcard and baits annotations for 'rgi bwt'.")

def load_reference_card_and_baits(local_db, card_fasta_file, baits_fasta_file, filename):
	db = get_location(local_db)

	filenames = []
	filenames.append(card_fasta_file)
	filenames.append(baits_fasta_file)

	# combine the two fastas
	import fileinput
	with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
		for line in fin:
			fout.write(line)
	# load baits only file
	load_reference_baits_only(local_db, baits_fasta_file)
	logger.info("loaded card and baits annotations for 'rgi bwt'.")

def load_reference_card_and_wildcard(local_db, card_fasta_file, wildcard_fasta_file, filename):

	db = get_location(local_db)
	# check files
	if card_fasta_file is None:
	# if os.path.isfile(os.path.abspath(card_fasta_file)) == False:
		logger.error("missing card reference file")
		exit()
	if os.path.isfile(os.path.abspath(wildcard_fasta_file)) == False:
		logger.error("missing wildcard reference file")
		exit()

	filenames = []
	filenames.append(card_fasta_file)
	filenames.append(wildcard_fasta_file)

	# combine the two fastas
	import fileinput
	with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
		for line in fin:
			fout.write(line)
	logger.info("loaded the following files: {} annotations for 'rgi bwt'.".format(",".join(filenames)))

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
	logger.info("local_db: {}, filepath: {}, filename: {},  validate_json: {}".format(
		local_db,
		filepath,
		filename,
		validate_json
	))
	db = get_location(local_db)

	try:
		# copy new file
		shutil.copyfile(filepath, os.path.join(db, filename))
		logger.info("file {} loaded ok".format(filename))
	except Exception as e:
		logger.warning("failed to copy json file: {}".format(e))

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
    parser.add_argument('-i', '--card_json', required=True, help='must be a card database json file')
    parser.add_argument('--card_annotation', required=False, help="annotated reference FASTA for homolog models only, created using rgi card_annotation")
    parser.add_argument('--card_annotation_all_models', required=False, help="annotated reference FASTA which includes all models created using rgi card_annotation")
    parser.add_argument('--wildcard_annotation', required=False, help="annotated reference FASTA for homolog models only, created using rgi wildcard_annotation")
    parser.add_argument('--wildcard_annotation_all_models', required=False, help="annotated reference FASTA which includes all models created using rgi wildcard_annotation")
    parser.add_argument('--wildcard_index', required=False, help="wildcard index file (index-for-model-sequences.txt)")
    parser.add_argument('--wildcard_version', required=False, help="specify variants version used")
    parser.add_argument('--baits_annotation', required=False, help="annotated reference FASTA")
    parser.add_argument('--baits_index', required=False, help="baits index file (baits-probes-with-sequence-info.txt)")
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
