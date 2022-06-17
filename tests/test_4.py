import pytest
import os, json
import subprocess as sp
from app.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# Run all tests with
# pytest test_4.py -v -rxs --color=auto --durations=0
# or
# pytest test_4.py -v -rxs --color=auto --durations=0 -k "protein"

@pytest.fixture
def rgi():
	return MainBase(api=True)

def run_rgi(rgi, read_one, read_two, aligner, threads, output_file):
	parser = rgi.bwt_args()
	rgi.bwt_run(parser.parse_args([
		'--read_one', read_one,
		'--read_two', read_two,
		'--aligner', aligner,
		'--threads', threads,
		'--output_file', output_file,
		'--clean',
		'--debug'
    ]))

def check_gene(gene):
    if gene == "tet(Q)":
        return True
    else:
        return False

def test_rgi_bwt(rgi):
    read_one = "10_R1.fastq.gz"
    read_two = "10_R2.fastq.gz"
    filename = "output"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi(rgi, os.path.join(working_directory,inputs,read_one), os.path.join(working_directory,inputs,read_two), "bowtie2", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True