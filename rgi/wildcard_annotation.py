import os, sys, json, csv, argparse, glob
from Bio import SeqIO
from rgi.settings import *

def main(args):
	working_directory = os.getcwd()
	logger.info(json.dumps(args.__dict__, indent=2))
	version = args.version
	files = glob.glob(os.path.join(args.input_directory,"*"))
	tab_file = os.path.join(working_directory, "wildcard_annotations_v{}.txt".format(version))
	fasta_file = os.path.join(working_directory, "wildcard_database_v{}.fasta".format(version))
	all_fasta = os.path.join(working_directory, "wildcard_all.fasta")
	fastas = {}
	annotations = {}

	for f in files:
		if "index" in f:
			index_file = f
			with open(os.path.join(index_file), 'r') as ifile:
				lines = csv.reader(ifile, delimiter='\t')
				for row in lines:
					if row[0] != 'prevalence_sequence_id':
						annotations[row[0]] = {
							"model_id": row[1],
							"aro_term": row[2],
							"aro_accession":  row[3],
							"detection_model": row[4],
							"percent_identity": row[9],
							"drug_class": row[13],
							"resistance_mechanism": row[12],
							"amr_gene_family": row[11]
						}
	prev_models = get_model(args.input_directory)

	for f in files:
		if "nucleotide_fasta_protein_homolog_model_variants.fasta" in f:
			for record in SeqIO.parse(f, 'fasta'):
				if record:
					desc = record.description.replace(" ", "_")
					arr = desc.split("|")
					model_id = prev_models[arr[0]]
					Prevalence_Sequence_ID = arr[0].split(":")[-1]
					header = "Prevalence_Sequence_ID:{prevalence_sequence_id}|ID:{model_id}|Name:{aro_term}|{aro_accession}".format(
								prevalence_sequence_id=Prevalence_Sequence_ID,
								model_id=annotations[Prevalence_Sequence_ID]["model_id"],
								aro_term=annotations[Prevalence_Sequence_ID]["aro_term"].replace(" ", "_"),
								aro_accession=annotations[Prevalence_Sequence_ID]["aro_accession"]
							)
					fastas[header] = str(record.seq)

	"""
	write fasta (FASTA format)
	"""
	with open(fasta_file, 'w') as fout:
		for i in fastas:
			fout.write(">{}\n".format(i))
			fout.write("{}\n".format(fastas[i]))

	print("Done writing wildcard_database_v{}.fasta".format(version))

def get_model(input_directory):
	prev_models = {}
	try:
		with open(os.path.join(input_directory, "index-for-model-sequences.txt"), 'r') as ifile:
			lines = csv.reader(ifile, delimiter='\t')
			for row in lines:
				if row[0] != 'prevalence_sequence_id':
					prev_models["Prevalence_Sequence_ID:{}".format(row[0])] = row[1]
	except Exception as e:
		raise e
	else:
		pass
	return prev_models

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi wildcard_annotation", description='Creates card annotations for RGI BWT from Variants or Wilcard data')
    parser.add_argument('-i', '--input_directory', dest="input_directory", required=True, help="input directory for wildcard")
    parser.add_argument('-v', '--version', dest="version", required=True, help="specify version downloaded for wildcard / variants")
    parser.add_argument('-j', '--card_json', dest="card_json", required=True, help="card.json file")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
