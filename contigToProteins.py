import os
import sys

def catProteins(afile):

	with open(afile, 'r') as f:
		data = f.readlines()

	with open('contig.fsa', 'w') as wf:

		startrecord = False

		for eachline in data:
			#eachline = eachline.strip()
			if eachline == '':
				pass
			elif eachline[0:18] == "Model information:":
				startrecord = False
			elif startrecord:
				if eachline[0] == '>':
					print>>wf, eachline.replace('\t>', '|').strip()
				else:
					print>>wf, eachline.strip()
			elif eachline[0:19] == "Predicted proteins:":
				startrecord = True


def main(argvfile):
	os.system("./mgm/gmhmmp -r -m ./mgm/MetaGeneMark_v1.mod -o draft -a " + argvfile)
	catProteins('./draft')
	os.remove('./draft')


if __name__ == '__main__':
	main(sys.argv[1])
