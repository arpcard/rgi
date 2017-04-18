import os
import sys
import filepaths
import time
import hashlib
import datetime

script_path = filepaths.determine_path()
working_directory = os.getcwd()

path = script_path

def catProteins(filename,afile):

	with open(afile, 'r') as f:
		data = f.readlines()
	f.close()

	with open(working_directory+'/'+filename+'.contigToORF.fsa', 'w') as wf:

		startrecord = False

		for eachline in data:
			if eachline == '':
				pass
			elif eachline[0:18] == "Model information:":
				startrecord = False
			elif startrecord:
				if eachline[0] == '>':
					print>>wf, eachline.replace('\t>', '|').strip()
				else:
					if eachline.strip() != "":
						print>>wf, eachline.strip()
			elif eachline.strip() == "Nucleotide sequence of predicted genes:":
				startrecord = True
	wf.close()

def get_character_len(file_name):
	chars = words = lines = 0
	with open(file_name, 'r') as in_file:
	    for line in in_file:
	    	if line[0] == '>':
	    		pass
	    	else:
		        lines += 1
		        words += len(line.split())
		        chars += len(line)
	in_file.close()
	return chars	

def main(argvfile,clean,orf):
	clean = clean.lower()
	orf = orf.lower()
	if orf == "genemark":
		orf_metagenemark(argvfile,clean)
	else:
		orf_prodigal(argvfile,clean)
   

def orf_metagenemark(argvfile,clean):
	filename = os.path.basename(argvfile)
	os.system(path+"/mgm/gmhmmp -r -m "+path+"/mgm/MetaGeneMark_v1.mod -o "+working_directory+"/"+filename+".adraft -d " + argvfile)
	catProteins(filename,working_directory +"/"+filename+".adraft")
	if clean == "yes":
		os.remove(working_directory +"/"+filename+".adraft")

'''
Training Using genomes:

Train:
	prodigal -i Salmonella_enterica_genome.fasta -p single -t Salmonella_enterica.trn

Using traning file:

	prodigal -i query.fasta -c -m -a trans_proteins.fasta -d trans_nulc -o output.txt -t Salmonella_enterica.trn -p single
'''
# prodigal -a dna.fasta.adraft -i dna.fasta -o dna.fasta.draft -p meta -d dna.fasta.contigToORF.fsa
def orf_prodigal(argvfile,clean):

	p = "single"
	count =  int(get_character_len(argvfile))
	#check for number of characters for selection procedure -- ribosomal binding sites
	if count < 200000:
		p = "meta"

	filename = os.path.basename(argvfile)
	os.system("prodigal -c -m -a "+working_directory+"/"+filename+".contig.fsa -i "+argvfile+" -o "+working_directory+"/"+filename+".draft -p "+p+" -d "+working_directory+"/"+filename+".contigToORF.fsa -q")

	# format the contig file headers to remove space
	format_fasta_headers(working_directory+"/"+filename+".contig.fsa")

	if clean == "yes":
		os.remove(working_directory +"/"+filename+".draft")

def getHashName(name):
	m = hashlib.md5()
	t = datetime.datetime.utcnow()
	m.update(name + str(t))
	return m.hexdigest()

def format_fasta_headers(afile):
	name = getHashName(afile)
	ofile = open(name, "w")
	from Bio import SeqIO
	for record in SeqIO.parse(afile, 'fasta'):
		new_header =  record.description.strip().replace(" # ", "#")
		ofile.write(">" + new_header + "\n" + str(record.seq) + "\n")
	ofile.close()
	os.rename(name, afile)	


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2], sys.argv[3])
