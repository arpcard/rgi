import os, sys, json, csv, argparse
"""
This script it used to create annotations and fasta for AMR++ using card data (version 2.0.0 and up)
"""

def main(args):
	working_directory = os.getcwd()
	"""
	reading card.json
	"""
	with open(os.path.join(args.input), 'r') as jfile:
		data = json.load(jfile)

	# get version
	try:
		version = data["_version"]
	except Exception as e:
		print("Error: missing version number")
		exit()

	annotations = []
	"""
	write card reference fasta (FASTA format)
	"""
	with open(os.path.join(working_directory, "card_database_v{}.fasta".format(version)), 'w') as fout:
		for i in data:
			if i.isdigit():
				# use homolog models only
				if data[i]['model_type_id'] in ['40292']:

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
							# print(data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"])

							if args.ncbi == True:
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
						print("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))
						
	"""
	write card anntotation fasta (tab-delimited format)
	"""
	# with open(os.path.join(working_directory, "card_annotations_v{}.txt".format(version)), "w") as af:
	# 	writer = csv.writer(af, delimiter=',', dialect='excel')
	# 	writer.writerow(["header","class", "mechanism", "group"])
	# 	for row in annotations:
	# 		writer.writerow(row) 	
		

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi card_annotation",description='Creates card annotations for RGI BWT from card.json')
    parser.add_argument('-i', '--input', dest="input", required=True, help="card.json file")
    parser.add_argument('--ncbi', dest="ncbi", action="store_true", help="adds ncbi accession to FASTA headers")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()

