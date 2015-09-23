import sys


def convertFqToFsa(managedfile):
	with open(managedfile, 'r') as f:
		data = f.readlines()

	linenum = 0
	'''	if linenum == 0: @description
	   	elif linenum == 1: DNA_Sequence
		elif linenum == 2: +More_description
		elif linenum == 3: Quality(Must be Same length as linenum[1])'''

	lineloc = 0

	for eachline in data:
		if linenum == 3:
			linenum = 0
		elif linenum == 2 and eachline[0] == '+':
			linenum = 3
		elif linenum == 1:
			print eachline
			linenum = 2
		elif eachline[0] == '@':
			linenum = 1
			print ('>' + eachline[1:].rstrip() + ' ' + data[lineloc+2][1:].rstrip())

		lineloc += 1


def main(argvfile):
	sys.stdout = open('read.fsa', 'w')
	convertFqToFsa(argvfile)


if __name__ == "__main__":
    main(sys.argv[1])
