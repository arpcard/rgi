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


def printCSV(resultfile,ofile,orf,verbose):
	if os.path.isfile(resultfile) == False:
		print>>sys.stderr, "convertJsonToTSV missing input JSON file."
		exit()

	try:
		with open(resultfile, 'r') as f:
			data = json.load(f)
		f.close()
	
	except ValueError:
		print>>sys.stderr, "convertJsonToTSV expects a file contains a VALID JSON string."
		exit()

	with open(working_directory+"/"+ofile+".txt", "w") as af:
		writer = csv.writer(af, delimiter='\t', dialect='excel')
		writer.writerow(["ORF_ID", "CONTIG", "START", "STOP", "ORIENTATION", "CUT_OFF", "PASS_EVALUE", "Best_Hit_evalue", "Best_Hit_ARO", "Best_Identities", "ARO", "ARO_name", "Model_type", "SNP", "Best_Hit_ARO_category", "ARO_category", "PASS_bitscore", "Best_Hit_bitscore", "bit_score","Predicted_DNA","Predicted_Protein","CARD_Protein_Sequence","LABEL","ID","Model_id"])
		for item in data:
			minevalue = 0.0
			minscore = 0.0
			maxpercent = 0.0
			startCompare = False
			minARO = ""
			bestAROcategorydict = {}
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
			predictedDNA = ""
			geneID = ""
			hitID = ""
			topModel = ""
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
					"""
					if data[item][it]["model_type_id"] == 41091:
						if checkKeyExisted("SNP",data[item][it]):
							temp = data[item][it]["SNP"]["original"] + str(data[item][it]["SNP"]["position"]) + data[item][it]["SNP"]["change"]
							snpList.append(convert(temp))
						else:
							snpList.append("n/a")
					"""

					AROlist.append(convert(data[item][it]["ARO_accession"]))
					AROnameList.append(convert(data[item][it]["ARO_name"]))
					bitScoreList.append(data[item][it]["bit-score"])
					pass_evalue = str(data[item][it]["pass_evalue"]).split("|")[0]
					pass_bitscore = "n/a"
					AROcatList.append(cgList)
					typeList.append(convert(data[item][it]["model_type"]))
					cutoffList.append(convert(data[item][it]["type_match"]))
					identityList.append(float(data[item][it]["perc_identity"]))

					bestAROcategory = []

					# sort results by minimum e-value and maximum percent identity
					if startCompare:
						if maxscore < data[item][it]["bit-score"] and maxpercent < float(data[item][it]["perc_identity"]):
							minevalue = data[item][it]["evalue"]
							maxscore = data[item][it]["bit-score"]
							maxpercent = float(data[item][it]["perc_identity"])
							minARO = data[item][it]["ARO_name"]
							topModel = data[item][it]["model_id"]
							SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]

							if "orf_prot_sequence" in data[item][it]:
								predictedProtein = data[item][it]["orf_prot_sequence"]
							if "orf_dna_sequence" in data[item][it]:
								predictedDNA = data[item][it]["orf_dna_sequence"]

							if checkKeyExisted("ARO_category", data[item][it]):
								for key in data[item][it]["ARO_category"]:
									bestAROcategory.append(str(data[item][it]["ARO_category"][key]["category_aro_name"].encode('ascii','replace')))
								bestAROcategorydict[str(minARO)+"|"+str(minevalue)] = bestAROcategory

							if "hsp_num:" in it:
								hitID = it
						

					else:
						startCompare = True
						minevalue = data[item][it]["evalue"]
						maxscore = data[item][it]["bit-score"]
						maxpercent = float(data[item][it]["perc_identity"])
						minARO = data[item][it]["ARO_name"]
						topModel = data[item][it]["model_id"]
						SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]

						if "orf_prot_sequence" in data[item][it]:
							predictedProtein = data[item][it]["orf_prot_sequence"]
						if "orf_dna_sequence" in data[item][it]:
								predictedDNA = data[item][it]["orf_dna_sequence"]

						if checkKeyExisted("ARO_category", data[item][it]):
							for key in data[item][it]["ARO_category"]:
								bestAROcategory.append(str(data[item][it]["ARO_category"][key]["category_aro_name"].encode('ascii','replace')))
							bestAROcategorydict[str(minARO)+"|"+str(minevalue)] = bestAROcategory

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
					if orf == "genemark":
						#for protein RGI runs where there's no | or seq_start/stop/strand
						if findnthbar(item, 4) == "":
							writer.writerow([item, 
								"", 
								"", 
								"", 
								"", 
								', '.join(list(clist)),
								pass_evalue,
								minevalue, 
								minARO, 
								maxpercent, 
								', '.join(map(lambda x:"ARO:"+x, AROlist)), 
								'; '.join(list(arocatset)),
								'; '.join(list(tl)), 
								snpList, 
								'; '.join(bestAROcategorydict[str(minARO)+"|"+str(minevalue)]) ,
								'; '.join(AROsortedList),
								pass_bitscore,
								maxscore ,
								', '.join(map(str, bitScoreList)),
								predictedDNA,
								predictedProtein,
								SequenceFromBroadStreet,
								geneID,
								hitID, 
								topModel
								])
		                                else:
						        writer.writerow([findnthbar(item, 0), 
						        	findORFfrom(item), 
						        	int(findnthbar(item, 4))-1, 
						        	int(findnthbar(item, 5))-1, 
						        	findnthbar(item, 3), 
						        	', '.join(list(clist)), 
						        	pass_evalue,
						        	minevalue ,
						        	minARO, 
						        	max(identityList), 
						        	', '.join(map(lambda x:"ARO:"+x, AROlist)), 
						        	'; '.join(list(arocatset)), 
						        	'; '.join(list(tl)), 
						        	snpList, 
						        	'; '.join(bestAROcategorydict[str(minARO)+"|"+str(minevalue)]) ,
						        	'; '.join(AROsortedList), 
						        	pass_bitscore,
						        	maxscore ,
						        	', '.join(map(str, bitScoreList)),
						        	predictedDNA,
						        	predictedProtein,
						        	SequenceFromBroadStreet,
						        	geneID,
						        	hitID,
						        	topModel
						        	])
					else:
						if findnthbar2(item, 1) == "":
							writer.writerow([item, 
								"", 
								"", 
								"", 
								"", 
								', '.join(list(clist)),
								pass_evalue,
								minevalue,
								minARO, 
								maxpercent, 
								', '.join(map(lambda x:"ARO:"+x, AROlist)), 
								'; '.join(list(arocatset)),
								', '.join(list(tl)), 
								snpList, 
								'; '.join(bestAROcategorydict[str(minARO)+"|"+str(minevalue)]),
								'; '.join(AROsortedList), 
								pass_bitscore,
								maxscore,
								', '.join(map(str, bitScoreList)),
								predictedDNA,
								predictedProtein,
								SequenceFromBroadStreet,
								geneID,
								hitID,
								topModel
								])
						else:
							writer.writerow([findnthbar2(item, 0),
								findnthbar2(item, 4).strip(" "), 
					        	int(findnthbar2(item, 1))-1, 
					        	int(findnthbar2(item, 2))-1, 
					        	findnthbar2(item, 3), 
					        	', '.join(list(clist)), 
					        	pass_evalue, 
					        	minevalue,
					        	minARO, 
					        	maxpercent, 
					        	', '.join(map(lambda x:"ARO:"+x, AROlist)), 
					        	', '.join(list(arocatset)), 
					        	', '.join(list(tl)), 
					        	snpList, 
					        	'; '.join(bestAROcategorydict[str(minARO)+"|"+str(minevalue)]),
					        	'; '.join(AROsortedList), 
					        	pass_bitscore,
					        	maxscore,
					        	', '.join(map(str, bitScoreList)),
					        	predictedDNA,
					        	predictedProtein,
					        	SequenceFromBroadStreet,
					        	geneID,
					        	hitID,
					        	topModel
						        ])



	af.close()

def manual():
	h = {}
	h["ORF_ID"] = "Open Reading Frame identifier (internal to RGI)"
	h["CONTIG"] = "Source Sequence"
	h["START"] = "Start co-ordinate of ORF"
	h["STOP"] = "End co-ordinate of ORF"
	h["ORIENTATION"] = "Strand of ORF"
	h["CUT_OFF"] = "RGI Detection Paradigm"
	h["PASS_EVALUE"] = "STRICT detection model Expectation value cut-off"
	h["Best_Hit_evalue"] = "Expectation value of match to top hit in CARD"
	h["Best_Hit_ARO"] = "ARO term of top hit in CARD"
	h["Best_Identities"] = "Percent identity of match to top hit in CARD"
	h["ARO"] = "ARO accession of top hit in CARD"
	h["ARO_name"] = "ARO term of top hit in CARD"
	h["Model_type"] = "CARD detection model type"
	h["SNP"] = "Observed mutation (if applicable)"
	h["Best_Hit_ARO_category"] = "top hit ARO Categorization"
	h["ARO_category"] = "ARO Categorization"
	h["PASS_bitscore"] = "STRICT detection model bitscore value cut-off"
	h["Best_Hit_bitscore"] = "Bit score of match to top hit in CARD"
	h["bit_score"] = "Bitscore of match to top hit in CARD"
	h["Predicted_DNA"] = "ORF predicted nucleotide sequence"
	h["Predicted_Protein"] = "ORF predicted protein sequence"
	h["CARD_Protein_Sequence"] = "Protein sequence of top hit in CARD"
	h["LABEL"] = "ORF label (internal to RGI)"
	h["ID"] = "HSP identifier (internal to RGI)"
	h["Model_id"] = "CARD detection model id"

	print "\n"
	print "COLUMN","\t\t\t","HELP_MESSAGE"
	for i in h:
		print i,"\t\t\t",h[i]	
	print "\n"


class customAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        manual()
        exit()

def main(args):
	afile = args.afile
	ofile = args.output
	#orf = args.orf.lower()
	orf = "prodigal"
	verbose = args.verbose.lower()

	# Check if file is compressed
	if afile.endswith('.gz'):
		afile = rgi.decompress(afile,'gz',working_directory)

	if os.path.isfile(afile):	
		printCSV(afile,ofile,orf,verbose)
	else:
		print "Missing file: ",afile 
	rgi.removeTemp()

def run():
	parser = argparse.ArgumentParser(description='Convert RGI JSON file to Tab-delimited file')
	parser.add_argument('-i','--afile',help='must be a json file generated from RGI in JSON or gzip format e.g out.json, out.json.gz')	
	parser.add_argument('-o', '--out_file',  dest="output", default="dataSummary", help="Output Tab-delimited file (default=dataSummary)")
	parser.add_argument('-v', '--verbose', dest="verbose", default="OFF", help = "include help menu. Options are OFF or ON  (default = OFF for no help)")
	parser.add_argument('--headers', dest="headers", action=customAction,nargs=0,  help = "print tab-delimted help. Options are OFF or ON  (default = OFF for no help)")
	args = parser.parse_args()
	main(args)	

if __name__ == '__main__':
	run()

