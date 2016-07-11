import json
import csv
import sys
import os
import rgi
import re

import argparse
import filepaths

script_path = filepaths.determine_path()
working_directory = os.getcwd()
path = script_path+"/"

def check_delimiter(fastaHeader):
	# Colon
	if ((':').join(re.split(':',fastaHeader)) == fastaHeader) and re.split(':',fastaHeader)[0] != fastaHeader:
		return ":"
	# Pipe
	elif (('|').join(re.split('|',fastaHeader)) == fastaHeader) and re.split('|',fastaHeader)[0] != fastaHeader:
		return "|"
	# Dash
	elif (('-').join(re.split('-',fastaHeader)) == fastaHeader) and re.split('-',fastaHeader)[0] != fastaHeader:
		return "-"
	# Underscore
	elif (('_').join(re.split('_',fastaHeader)) == fastaHeader) and re.split('_',fastaHeader)[0] != fastaHeader:
		return "_"	
	# Other
	else:
		return ""	

#output the information particular field from alignment.Title by splicing it by '|'
def findnthbar(bunchstr, start):

	barc = 0
	over = start+1
	temp = ""

	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if barc == start:
			if eachc == '|':
				pass
			else:
				temp += eachc
		if barc == over:
			break		

	return temp

#output the information particular field from alignment.Title by splicing it by '#'
def findnthbar2(bunchstr, n):
	arr = bunchstr.split("#")
	if n < len(arr):
		# gene id
		if n == 1 and arr[n]:
			return int(arr[n])
		elif n == 2:
			return int(arr[n])
		elif n == 3:
			 if int(arr[n]) == 1:
			 	# positive
			 	return "+"
			 else:
			 	# neg
			 	return "-"
		else:
			return arr[n]
	else:
		return ""

def findORFfrom (bunchstr):
	barc = 0
	start = 6
	temp = ""
	allout = False

	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if allout or barc == start:
			allout = True
			temp += eachc

	return temp[1:]


def convert(input):
	if isinstance(input, dict):
		return dict((convert(key), convert(value)) for key, value in input.iteritems())
	elif isinstance(input, list):
		return [convert(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


def printCSV(resultfile,ofile,orf):
	try:
		with open(resultfile, 'r') as f:
			data = json.load(f)
	
	except ValueError:
		print>>sys.stderr, "convertJsonToTSV expects a file contains a VALID JSON string."
		exit()

	with open(working_directory+"/"+ofile+".txt", "w") as af:
		writer = csv.writer(af, delimiter='\t', dialect='excel')
		writer.writerow(["ORF_ID", "CONTIG", "START", "STOP", "ORIENTATION", "CUT_OFF", "PASS_EVALUE", "Best_Hit_evalue", "Best_Hit_ARO", "Best_Identities", "ARO", "ARO_name", "Model_type", "SNP", "AR0_category", "bit_score","Predicted_Protein","CARD_Protein_Sequence","LABEL","ID"])
		for item in data:
			minevalue = False
			startCompare = False
			minARO = 0
			AROlist = []
			AROnameList = []
			bitScoreList = []
			AROcatList = []
			snpList = []
			cutoffList = []
			typeList = []
			evalueList = []
			identityList = []
			SequenceFromBroadStreet = ""
			predictedProtein = ""
			geneID = ""
			hitID = ""

			if item not in ["_metadata","data_type"]:
				geneID = item

				for it in data[item]:
					cgList = []
					if checkKeyExisted("ARO_category", data[item][it]):
						for aroctkey in data[item][it]["ARO_category"]:
							cgList.append(str(data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace')))

					if data[item][it]["model_type_id"] == 40293:
						temp = data[item][it]["SNP"]["original"] + str(data[item][it]["SNP"]["position"]) + data[item][it]["SNP"]["change"]
						snpList.append(convert(temp))
					elif data[item][it]["model_type_id"] == 40292:
						snpList.append("n/a")

					AROlist.append(convert(data[item][it]["ARO_accession"]))
					AROnameList.append(convert(data[item][it]["ARO_name"]))
					bitScoreList.append(data[item][it]["bit-score"])
					pass_evalue = str(data[item][it]["pass_evalue"]).split("|")[0]
					AROcatList.append(cgList)
					typeList.append(convert(data[item][it]["model_type"]))
					cutoffList.append(convert(data[item][it]["type_match"]))
					idenPercent = float(data[item][it]["max-identities"]) / len(data[item][it]["query"])
					'''print>>sys.stderr, data[item][it]["max-identities"]
					print>>sys.stderr, len(data[item][it]["query"])
					print (str(269/289) + "haha")
					print>>sys.stderr, float(data[item][it]["max-identities"] % len(data[item][it]["query"]))'''
					identityList.append(idenPercent)
					
					if startCompare:
						if minevalue > data[item][it]["evalue"]:
							minevalue = data[item][it]["evalue"]
							minARO = data[item][it]["ARO_name"]
							SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]
							if "orf_prot_sequence" in data[item][it]:
								predictedProtein = data[item][it]["orf_prot_sequence"]
							if "hsp_num:" in it:
								hitID = it							
					else:
						startCompare = True
						minevalue = data[item][it]["evalue"]
						minARO = data[item][it]["ARO_name"]
						SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]
						if "orf_prot_sequence" in data[item][it]:
							predictedProtein = data[item][it]["orf_prot_sequence"]
						if "hsp_num:" in it:
							hitID = it						

			clist = set(cutoffList)
			tl = set(typeList)
			arocatset = set(AROnameList)
			
			if set(snpList) == set(['n/a']):
				snpList = 'n/a'
			else:
				snpList = ', '.join(snpList)

			from itertools import chain
			AROcatList = list(chain.from_iterable(AROcatList))
			AROcatalphaSet = set(AROcatList)
			AROsortedList = sorted(list(AROcatalphaSet))

			if typeList:
				if orf == "1":
					#for protein RGI runs where there's no | or seq_start/stop/strand
					if findnthbar(item, 4) == "":
						writer.writerow([item, "", "", "", "", ', '.join(list(clist)),pass_evalue, minevalue, minARO, max(identityList), ', '.join(map(lambda x:"ARO:"+x, AROlist)), ', '.join(list(arocatset)),', '.join(list(tl)), snpList, ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)),predictedProtein,SequenceFromBroadStreet,geneID,hitID])
	                                else:
					        writer.writerow([findnthbar(item, 0), findORFfrom(item), int(findnthbar(item, 4))-1, int(findnthbar(item, 5))-1, findnthbar(item, 3), ', '.join(list(clist)), pass_evalue, minevalue, minARO, max(identityList), ', '.join(map(lambda x:"ARO:"+x, AROlist)), ', '.join(list(arocatset)), ', '.join(list(tl)), snpList, ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)),predictedProtein,SequenceFromBroadStreet,geneID,hitID])
				else:
					if findnthbar2(item, 1) == "":
						writer.writerow([item, "", "", "", "", ', '.join(list(clist)),pass_evalue, minevalue, minARO, max(identityList), ', '.join(map(lambda x:"ARO:"+x, AROlist)), ', '.join(list(arocatset)),', '.join(list(tl)), snpList, ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)),predictedProtein,SequenceFromBroadStreet,geneID,hitID])
	                                else:
					        writer.writerow([findnthbar2(item, 0), 
					        	findnthbar2(item, 4).strip(" "), 
					        	int(findnthbar2(item, 1))-1, 
					        	int(findnthbar2(item, 2))-1, 
					        	findnthbar2(item, 3), 
					        	', '.join(list(clist)), pass_evalue, minevalue, minARO, max(identityList), ', '.join(map(lambda x:"ARO:"+x, AROlist)), ', '.join(list(arocatset)), ', '.join(list(tl)), snpList, ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)),predictedProtein,SequenceFromBroadStreet,geneID,hitID])



def main(args):
	afile = args.afile
	ofile = args.output
	orf = args.orf
	# Check if file is compressed
	if afile.endswith('.gz'):
		afile = rgi.decompress(afile,'gz',working_directory)

	printCSV(afile,ofile,orf)
	rgi.removeTemp()

def run():
	"""required: sys.argv[1] must be a json file"""
	parser = argparse.ArgumentParser(description='Convert RGI JSON file to Tab-delimited file')
	parser.add_argument('-i','--afile',help='must be a json file generated from RGI in JSON or gzip format e.g out.json, out.json.gz')	
	parser.add_argument('-o', '--out_file',  dest="output", default="dataSummary", help="Output JSON file (default=dataSummary)")
	parser.add_argument('-x', '--orf', dest="orf", default="0", help = "choose between prodigal and MetaGeneMark orf finder. Options are 0 or 1  (default = 0 for using prodigal)")
	args = parser.parse_args()
	main(args)	

if __name__ == '__main__':
	run()

