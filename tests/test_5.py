import pytest
import os, json
import subprocess as sp
from app.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# Run all tests with
# pytest test_5.py -v -rxs --color=auto --durations=0
# or
# pytest test_.py -v -rxs --color=auto --durations=0 -k "protein"

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

def run_rgi_single(rgi, read_one, aligner, threads, output_file):
    parser = rgi.bwt_args()
    rgi.bwt_run(parser.parse_args([
        '--read_one', read_one,
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

def test_rgi_bwt_bwa(rgi):
    read_one = "10_R1.fastq.gz"
    read_two = "10_R2.fastq.gz"
    filename = "output_bwt"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi(rgi, os.path.join(working_directory,inputs,read_one), os.path.join(working_directory,inputs,read_two), "bwa", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True

def test_rgi_bwt_bwa_single(rgi):
    read_one = "10_R1.fastq.gz"
    filename = "output_bwt_bwa_single"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi_single(rgi, os.path.join(working_directory,inputs,read_one), "bwa", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True

def test_rgi_bwt_bowtie2(rgi):
    read_one = "10_R1.fastq.gz"
    read_two = "10_R2.fastq.gz"
    filename = "output_bwt_bowtie2"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi(rgi, os.path.join(working_directory,inputs,read_one), os.path.join(working_directory,inputs,read_two), "bowtie2", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True

def test_rgi_bwt_bowtie2_single(rgi):
    read_one = "10_R1.fastq.gz"
    filename = "output_bwt_bowtie2_single"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi_single(rgi, os.path.join(working_directory,inputs,read_one), "bowtie2", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True

def test_rgi_bwt_kma(rgi):
    read_one = "10_R1.fastq.gz"
    read_two = "10_R2.fastq.gz"
    filename = "output_bwt_kma"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi(rgi, os.path.join(working_directory,inputs,read_one), os.path.join(working_directory,inputs,read_two), "kma", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True


def test_rgi_bwt_kma_interleaved(rgi):
    read_one = "10_R1.fastq.gz"
    filename = "output_bwt_kma_interleaved"
    output_file = os.path.join(working_directory,outputs,"{}".format(filename))
    run_rgi_single(rgi, os.path.join(working_directory,inputs,read_one), "kma", "3" ,output_file)
    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{amr_gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        amr_gene = "tet(Q)"
    ))
    assert check_gene(gene) == True

    
