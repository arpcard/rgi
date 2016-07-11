import os
import sys
import filepaths

script_path = filepaths.determine_path()
working_directory = os.getcwd()

path = script_path

def catProteins(filename,afile):

	with open(afile, 'r') as f:
		data = f.readlines()

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
					print>>wf, eachline.strip()
			elif eachline.strip() == "Nucleotide sequence of predicted genes:":
				startrecord = True

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
	return chars	

def main(argvfile,clean,orf):
	if orf == "1":
		orf_metagenemark(argvfile,clean)
	else:
		orf_prodigal(argvfile,clean)
   

def orf_metagenemark(argvfile,clean):
	filename = os.path.basename(argvfile)
	os.system(path+"/mgm/gmhmmp -r -m "+path+"/mgm/MetaGeneMark_v1.mod -o "+working_directory+"/"+filename+".adraft -d " + argvfile)
	catProteins(filename,working_directory +"/"+filename+".adraft")
	if clean == "1":
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

	#check for number of characters for selection procedure -- ribosomal binding sites
	if get_character_len(argvfile) < 20000:
		p = "meta"

	filename = os.path.basename(argvfile)
	os.system("prodigal -c -m -a "+working_directory+"/"+filename+".contig.fsa -i "+argvfile+" -o "+working_directory+"/"+filename+".draft -p "+p+" -d "+working_directory+"/"+filename+".contigToORF.fsa -q")

	if clean == "1":
		os.remove(working_directory +"/"+filename+".draft")

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2], sys.argv[3])
