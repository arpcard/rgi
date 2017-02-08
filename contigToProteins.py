import os
import sys
import filepaths

script_path = filepaths.determine_path()
working_directory = os.getcwd()

path = script_path

def catProteins(filename,afile):

	with open(afile, 'r') as f:
		data = f.readlines()
	f.close()
	
	with open(working_directory+'/'+filename+'.contig.fsa', 'w') as wf:

		startrecord = False

		for eachline in data:
			if eachline == '':
				pass
			elif eachline[0:18] == "Model information:":
				startrecord = False
			elif startrecord:
				if eachline[0] == '>':
					print>>wf, eachline.replace('\t>', '|').strip().lstrip()
				else:
					if eachline.strip() != "":
						print>>wf, eachline.strip()

			elif eachline[0:19] == "Predicted proteins:":
				startrecord = True
	wf.close()


def main(argvfile, clean):
	filename = os.path.basename(argvfile)
	os.system(path+"/mgm/gmhmmp -r -m "+path+"/mgm/MetaGeneMark_v1.mod -o "+working_directory+"/"+filename+".draft -a " + argvfile)
	catProteins(filename,working_directory +"/"+filename+".draft")
	if clean == "yes":
		os.remove(working_directory +"/"+filename+".draft")


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
