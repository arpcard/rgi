import pytest
import os, json
from app.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# Run all tests with
# pytest test_1.py -v -rxs --color=auto --durations=0
# or
# pytest test_1.py -v -rxs --color=auto --durations=0 -k "protein"

@pytest.fixture
def rgi():
	return MainBase(api=True)
"""
def import_from_ncbi(accession, seq_type="nucleotide"):
	from Bio import Entrez
	import time
	import json
	Entrez.email ="johnQ@gmail.com"
	try:
		handle = Entrez.efetch(db=seq_type, id=accession, retmode="xml" , rettype = "fasta")
	except Exception as e:
		time.sleep(20)
		handle = Entrez.efetch(db=seq_type, id=accession, retmode="xml" , rettype = "fasta")
	records = Entrez.read(handle)

	return records[0]["TSeq_sequence"]
"""

def validate_results(filepath,perc_identity=0, ARO_name='', type_match=''):
	pi = ""
	name = ""
	tm = ""
	filename = os.path.basename(filepath)
	f = os.path.join("{}".format(filepath))
	if os.path.isfile(f):
		with open(f) as json_file:
			json_data = json.load(json_file)
			for i in json_data:
				if i not in ["_metadata"]:
					for j in json_data[i]:
						for k in json_data[i][j]:
							pi = json_data[str(i)][str(j)]["perc_identity"]
							name = json_data[str(i)][str(j)]["ARO_name"]
							tm = json_data[str(i)][str(j)]["type_match"]
						if pi == perc_identity and name == ARO_name and tm == type_match:
							#print(pi, name, tm)
							return True
			return False
	else:
		print("missing file: {}".format(f))
		return False

def run_rgi(rgi, input_type, input_sequence, output_file):
	parser = rgi.main_args()
	rgi.main_run(parser.parse_args([
		'--input_type', input_type,
		'--input_sequence', input_sequence,
		'--output_file', output_file,
		'--alignment_tool', alignment_tool,
		'--clean',
		'--include_loose',
		'--include_nudge',
		'--low_quality',
		'--debug'
    ]))

# def run_rgi_with_nudge(rgi, input_type, input_sequence, output_file):
# 	parser = rgi.main_args()
# 	rgi.main_run(parser.parse_args([
# 		'--input_type', input_type,
# 		'--input_sequence', input_sequence,
# 		'--output_file', output_file,
# 		'--alignment_tool', alignment_tool,
# 		'--clean',
# 		'--include_loose',
# 		'--include_nudge',
# 		'--low_quality',
# 		'--debug'
#     ]))

def test_rgi_protein_sequence(rgi):

	filename = "protein.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'protein', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 100, 'NDM-1', 'Perfect') == True

def test_rgi_nucleotide_sequence(rgi):

	filename = "NC_020818.1.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'contig', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 100, 'NDM-1', 'Perfect') == True
	assert validate_results(output_file, 100, "APH(3')-VI", 'Perfect') == True
	assert validate_results(output_file, 100, 'mphE', 'Perfect') == True
	assert validate_results(output_file, 100, 'msrE', 'Perfect') == True

def test_rgi_homolog_model(rgi):

	filename = "homolog.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'contig', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 100, 'NDM-1', 'Perfect') == True

def test_rgi_variant_model(rgi):

	filename = "variant.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'protein', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 99.88, 'Escherichia coli gyrB conferring resistance to aminocoumarin', 'Strict') == True

def test_rgi_overexpression_model(rgi):

	filename = "overexpression.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'protein', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 100, 'nalC', 'Strict') == True
	assert validate_results(output_file, 99.53, 'nalC', 'Strict') == True

def test_rgi_effluxpump_model(rgi):

	filename = "effluxpump.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'protein', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 100, 'MexA', 'Perfect') == True
	assert validate_results(output_file, 100, 'MexB', 'Perfect') == True
	assert validate_results(output_file, 100, 'OprM', 'Perfect') == True
	assert validate_results(output_file, 100, 'Pseudomonas aeruginosa CpxR', 'Perfect') == True
	assert validate_results(output_file, 99.82, 'MexB', 'Strict') == True
	assert validate_results(output_file, 99.59, 'OprM', 'Strict') == True
	assert validate_results(output_file, 98.12, 'nalC', 'Strict') == True
	assert validate_results(output_file, 97.28, 'MexR', 'Strict') == True
	assert validate_results(output_file, 89.12, 'MexR', 'Loose') == True

def test_rgi_rrna_model(rgi):

	filename = "rrna.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'contig', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 99.97, \
		'Streptococcus pneumoniae 23S rRNA mutation conferring resistance to macrolides and streptogramins antibiotics', 'Strict') == True


def test_rgi_nudge_loose_to_strict(rgi):

	filename = "loose_to_strict.fasta"
	output_file = os.path.join(working_directory,outputs,"{}.json".format(filename))
	run_rgi(rgi, 'contig', os.path.join(working_directory,inputs,filename), output_file)

	assert validate_results(output_file, 98.06, 'Escherichia coli EF-Tu mutants conferring resistance to Pulvomycin', 'Strict') == True

