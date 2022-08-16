import os, sys, json, csv, argparse
from app.settings import *
from argparse import RawTextHelpFormatter
from app.settings import APP_NAME, SOFTWARE_VERSION
"""
This script it used to create annotations and fasta for AMR++ using card data (version 2.0.0 and up)
"""
def main(args):
	if args.debug:
		logger.setLevel(10)
	"""
	reading card.json
	"""
	with open(os.path.join(args.input), 'r') as jfile:
		data = json.load(jfile)

	# get version
	try:
		version = data["_version"]
	except Exception as e:
		logger.error("missing version number")
		exit()

	# write FASTA with only homolog models
	write_fasta_annotation_file(data, version, args.ncbi)
	# write FASTA with homolog, variant, rRNA gene variant, overexpression, knockout models
	write_fasta_annotation_file(data, version, args.ncbi, True)

def write_fasta_annotation_file(data, version, ncbi, all_model_type=False):
	working_directory = os.getcwd()
	annotations = []
	# use homolog only
	selected_model_types = ['40292']
	output_file = os.path.join(working_directory, "card_database_v{}.fasta".format(version))
	if all_model_type is True:
		# use homolog, variant, rRNA gene variant, overexpression, knockout models
		selected_model_types = ['40292','40293','40295','41091','40354']
		output_file = os.path.join(working_directory, "card_database_v{}_all.fasta".format(version))

	"""
	write card reference fasta (FASTA format)
	"""
	with open(output_file, 'w') as fout:
		for i in data:
			if i.isdigit():
				# check for model types
				if data[i]['model_type_id'] in selected_model_types:

					drug_class = []
					mechanism = []
					group = []

					if "ARO_category" in data[i]:
						for c in data[i]["ARO_category"]:
							if "category_aro_class_name" in data[i]["ARO_category"][c]:
								if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
									drug_class.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
								if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
									mechanism.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
								if data[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
									group.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
					try:
						for seq in data[i]['model_sequences']['sequence']:
							if ncbi is True:
								# header used to be able to validate CARD sequences with genbank sequences
								header = ("gb|{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}".format(
									ARO_accession=data[i]['ARO_accession'],
									model_id=data[i]['model_id'],
									ARO_name=(data[i]['ARO_name']).replace(" ", "_"),
									ncbi=data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
									))
							else:
								header = ("ARO:{}|ID:{}|Name:{}|NCBI:{}".format(
									data[i]['ARO_accession'],
									data[i]['model_id'],
									(data[i]['ARO_name']).replace(" ", "_"),
									data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
									))

							sequence = data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["sequence"]
							fout.write(">{}\n".format(header))
							fout.write("{}\n".format(sequence))
							annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])

					except Exception as e:
						logger.debug(e)
						logger.warning("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))
	logger.info("Done writing {}".format(os.path.basename(output_file)))

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi card_annotation",description='{} - {} - CARD Annotation \n\nCreates card annotations for RGI BWT from card.json'.format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--input', dest="input", required=True, help="card.json file")
	parser.add_argument('--ncbi', dest="ncbi", action="store_true", help="adds ncbi accession to FASTA headers (default: False)")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode (default: False)")
	return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
