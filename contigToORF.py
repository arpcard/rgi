import os
import sys
import filepaths

script_path = filepaths.determine_path()
working_directory = os.getcwd()

path = script_path

def catProteins(filename,afile):

	with open(afile, 'r') as f:
		data = f.readlines()

	with open(working_directory+'/'+filename+'.contigtoORF.fsa', 'w') as wf:

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


def main(argvfile):
	filename = os.path.basename(argvfile)
	os.system(path+"/mgm/gmhmmp -r -m "+path+"/mgm/MetaGeneMark_v1.mod -o "+working_directory+"/"+filename+".adraft -d " + argvfile)
	catProteins(filename,working_directory +"/"+filename+".adraft")
	#os.remove(working_directory +"/"+filename+".adraft")

if __name__ == '__main__':
	main(sys.argv[1])
