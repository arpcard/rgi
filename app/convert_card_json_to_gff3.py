import sys
import os
import json
import argparse
import csv

"""
This script is used to convert CARD json file to GFF3 format.
The card.json file can be downloaded at https://card.mcmaster.ca/download under CARD data.

Usage:
	python convert_card_json_to_gff3.py -i card.json

"""
def main(args):
	if args.input_file == None:
		exit("Missing input card json file")

	with open(os.path.join(args.input_file), 'r') as jfile:
		j = json.load(jfile)
		with open(os.path.join("card.gff3"), "w") as af:
			sequences = []
			headers = []
			body = []
			writer = csv.writer(af, delimiter='\t')
			headers.append(['##gff-version 3.2.1'])
			for i in j:
				if i.isdigit():
					if "model_sequences" in j[i].keys():
						for k in j[i]["model_sequences"]["sequence"]:
							_source = "{}".format("CARD")	
							_type = "{}".format("CDS")
							_phase = "{}".format(".")
							_score = "{}".format(".")
							_seqid = "{gi}{aro}".format(gi="gi|"+j[i]["model_sequences"]["sequence"][k]["dna_sequence"]["accession"], aro="|ARO:" + j[i]["ARO_accession"])
							_start = "{}".format( int(j[i]["model_sequences"]["sequence"][k]["dna_sequence"]["fmin"]) + 1)
							_end =  "{}".format(j[i]["model_sequences"]["sequence"][k]["dna_sequence"]["fmax"])
							_strand = "{}".format(j[i]["model_sequences"]["sequence"][k]["dna_sequence"]["strand"])
							_attributes = "Name={name};Alias={aro}".format(
								name=j[i]["ARO_name"],
								aro="ARO:" + j[i]["ARO_accession"])

							headers.append(['##sequence-region '+_seqid+' '+_start+' '+_end])
							body.append([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes])
							sequences.append(format_fasta(str(_seqid), "{}".format(j[i]["model_sequences"]["sequence"][k]["dna_sequence"]["sequence"])))
			# headers
			for head_item in headers:
				af.write(head_item[0]+ '\n')
			# body
			for body_item in body:
				writer.writerow(body_item)
			# footer
			writer.writerow(["##FASTA"])
			for sequence in sequences:
				af.write(sequence)


def format_fasta(name, sequence):
    fasta_string = '>' + name + '\n' + sequence + '\n'
    return fasta_string

def run():
	parser = argparse.ArgumentParser(description='convert card json to gff3')
	parser.add_argument('-i','--input_file', dest="input_file", default=None, required=True, help='card.json input file')
	args = parser.parse_args()
	main(args)	

if __name__ == '__main__':
	run()