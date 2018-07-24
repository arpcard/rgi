import os, sys, json, csv, argparse, glob
from Bio import SeqIO
from Bio.SeqUtils import GC
from app.settings import *

def get_gc_length_genes(card_json):
	info = {}
	working_directory = os.getcwd()
	"""
	reading card.json
	"""
	with open(os.path.join(card_json), 'r') as jfile:
		data = json.load(jfile)

	annotations = []

	# with open(os.path.join(working_directory, "card_database_v{}.fasta".format(version)), 'w') as fout:
	for i in data:
		if i.isdigit():
			# use homolog and variants models only
			if data[i]['model_type_id'] in ['40292', '40293']:
				# print("what", data[i])
				# exit()
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
						header = ("ARO:{}|ID:{}|Name:{}".format(
							data[i]['ARO_accession'],
							data[i]['model_id'],
							(data[i]['ARO_name']).replace(" ", "_")
							))
						sequence = data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["sequence"]
						# fout.write(">{}\n".format(header))
						# fout.write("{}\n".format(sequence))
						# annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])
						# print(i, data[i])
						info[i] = {
							"gc": GC(sequence),
							"length": len(sequence),
							"nruns": sequence.count('NNN'),
							"sequence": sequence
							}

				except Exception as e:
					pass
					# print("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))

	return info

def main(args):
	working_directory = os.getcwd()
	info = get_gc_length_genes(args.card_json)
	logger.info(json.dumps(args.__dict__, indent=2))
	# exit("Stop")
	version = args.version
	files = glob.glob(os.path.join(args.input_directory,"*"))
	tab_file = os.path.join(working_directory, "wildcard_annotations_v{}.txt".format(version))
	fasta_file = os.path.join(working_directory, "wildcard_database_v{}.fasta".format(version))
	all_fasta = os.path.join(working_directory, "wildcard_all.fasta")
	fastas = {}

	prev_models = get_model(args.input_directory)
	# print(prev_models)
	# exit("Done?")
	'''
	with open(all_fasta, 'w') as fout:
		for f in files:
			if "index" in f:
				index_file = f
			elif "nucleotide_fasta" in f:
				for record in SeqIO.parse(f, 'fasta'):
					if record:
						desc = record.description.replace(" ", "_")
						arr = desc.split("|")
						# print(arr[1].split(":")[1])
						model_id = prev_models[arr[0]]						
						fastas["{}|{}|{}".format(arr[0],arr[1],arr[2])] = record.seq
						fout.write(">{}|ID:{}|Name:{}\n".format(
							arr[0],
							model_id,
							(arr[1].split(":")[1]).replace(" ", "_")
							))
						fout.write("{}\n".format(record.seq))
	
	'''
	for f in files:
		if "index" in f:
			index_file = f
		elif "nucleotide_fasta" in f:
			for record in SeqIO.parse(f, 'fasta'):
				if record:
					desc = record.description.replace(" ", "_")
					arr = desc.split("|")
					model_id = prev_models[arr[0]]						
					fastas["{}|{}|{}".format(arr[0],arr[1],arr[2])] = str(record.seq)

	annotations = []

	"""
	write fasta (FASTA format)
	"""
	new_fasta = {}
	# new_tsv = []
	outliers = []

	print("1/3 steps - Done writing fasta")
	try:
		with open(os.path.join(index_file), 'r') as ifile:
			lines = csv.reader(ifile, delimiter='\t')
			for row in lines:
				if row[0] != 'prevalence_sequence_id':
					anno, selector = get_headers(row[0], row[1], row[2], row[3], row[4], row[9], row[13], row[12], row[11])
					# if anno not in annotations:
						# annotations.append(anno)
						# append_to_tsv(tab_file, anno)
						# append_to_fasta(fasta_file, anno[0], fastas[selector])
						# new_tsv.append(anno)
					new_fasta[anno[0]] = fastas[selector]
					sequence = fastas[selector]
					gc = GC(sequence)
					nruns = sequence.count('NNN')
					seq_len = len(sequence)
					# sequence.count('NNN') > 0 or 
					try:
						if ((int(info[row[1]]["gc"]) - 20) <= int(gc) <= (int(info[row[1]]["gc"]) + 20)) == False:
							# print("anno =>", anno[0], "selector =>", selector)
							# print(info[row[1]])
							outliers.append({
									"model_id": row[1],
									"variant": selector,
									"reference_gc": info[row[1]]["gc"],
									"reference_length": info[row[1]]["length"],
									"variant_gc": gc,
									"variant_length": seq_len,
									"nruns": nruns
									# "sequence": sequence
									})
					except Exception as e:
						print(e)

						# exit()
			print("2/3 steps - Done annotations")
	except Exception as e:
		raise e
	else:
		pass

	with open("outliers.json", 'w') as ptr:
		json.dump(outliers, ptr)

	with open(fasta_file, 'w') as fout:
		for i in new_fasta:
			if "N" in new_fasta[i]:
				print("undetermined base skipping >{}".format(i))
			else:
				fout.write(">{}\n".format(i))
				fout.write("{}\n".format(new_fasta[i]))

	"""
	write card anntotation fasta (tab-delimited format)
	"""
	# with open(tab_file, "w") as af:
	# 	writer = csv.writer(af, delimiter=',', dialect='excel')
	# 	writer.writerow(["header","class", "mechanism", "group"])
	# 	added = []
	# 	for r in new_tsv:
	# 		writer.writerow(r)

	print("3/3 steps - Done writing wildcard_database_v{}.fasta".format(version))

def get_model(input_directory):
	prev_models = {}
	try:
		with open(os.path.join(input_directory, "index-for-model-sequences.txt"), 'r') as ifile:
			lines = csv.reader(ifile, delimiter='\t')
			for row in lines:
				# print(row[0], "=>", row[1])
				if row[0] != 'prevalence_sequence_id':
					prev_models["Prevalence_Sequence_ID:{}".format(row[0])] = row[1]
	except Exception as e:
		raise e
	else:
		pass
	return prev_models

def append_to_tsv(filepath, row):
	with open(filepath, "a") as af:
		writer = csv.writer(af, delimiter=',', dialect='excel')
		writer.writerow(row)	

def append_to_fasta(filepath, header, sequence):
	with open(filepath, 'a') as fout:
		if sequence:
			fout.write(">{}\n".format(header))
			fout.write("{}\n".format(sequence))

def get_sequence(selector, f):
	sequence = ""
	with open(f, "r") as handle:
		for record in SeqIO.parse(handle, 'fasta'):
			if record:
				if selector in record.description:
					sequence = record.seq
		return sequence


def get_headers(prevalence_sequence_id, model_id, aro_term, aro_accession, detection_model, percent_identity, drug_class, resistance_mechanism, amr_gene_family):
	# ARO:3003964|hp1181|protein homolog model|1|2430
	# header = "{aro_accession}|{aro_term}|{detection_model}|{prevalence_sequence_id}|{model_id}".format(
	# 	prevalence_sequence_id=prevalence_sequence_id, 
	# 	model_id=model_id, 
	# 	aro_term=aro_term.replace(" ", "_"),
	# 	aro_accession=aro_accession, 
	# 	detection_model=detection_model.replace(" ", "_")
	# )
	# ARO:3003964|hp1181|protein homolog model|1|2430
	# header = "{aro_accession}|{aro_term}|{detection_model}|Prevalence_Sequence_ID:{prevalence_sequence_id}|ID:{model_id}".format(
	# 	prevalence_sequence_id=prevalence_sequence_id, 
	# 	model_id=model_id, 
	# 	aro_term=aro_term.replace(" ", "_"),
	# 	aro_accession=aro_accession, 
	# 	detection_model=detection_model.replace(" ", "_")
	# )
	# Prevalence_Sequence_ID:1|ARO_Name:hp1181|ARO:3003964
	header = "Prevalence_Sequence_ID:{prevalence_sequence_id}|ID:{model_id}|Name:{aro_term}|{aro_accession}".format(
		prevalence_sequence_id=prevalence_sequence_id,
		model_id=model_id,
		aro_term=aro_term.replace(" ", "_"),
		aro_accession=aro_accession
	)
	# ARO:3003964|ID:2430|Name:hp1181|Prevalence_Sequence_ID:1
	# header = "{aro_accession}|ID:{model_id}|Name:{aro_term}|Prevalence_Sequence_ID:{prevalence_sequence_id}".format(
	# 	prevalence_sequence_id=prevalence_sequence_id, 
	# 	model_id=model_id, 
	# 	aro_term=aro_term.replace(" ", "_"),
	# 	aro_accession=aro_accession, 
	# 	detection_model=detection_model.replace(" ", "_")
	# )
	# Prevalence_Sequence_ID:1|ARO_Name:hp1181|ARO:3003964
	selector = "Prevalence_Sequence_ID:{prevalence_sequence_id}|ARO_Name:{aro_term}|{aro_accession}".format(
		prevalence_sequence_id=prevalence_sequence_id,
		aro_term=aro_term.replace(" ", "_"),
		aro_accession=aro_accession
	)

	return [header, drug_class, resistance_mechanism, amr_gene_family], selector

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