import json
import csv
import sys
#import tryrgi


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


def printCSV(resultfile):
	try:
		with open(resultfile, 'r') as f:
			data = json.load(f)
	
	except ValueError:
		print>>sys.stderr, "convertJsonToTSV expects a file contains a VALID JSON string."
		exit()

	with open("dataSummary.txt", "w") as af:
		writer = csv.writer(af, delimiter='\t', dialect='excel')
		writer.writerow(["ORF_ID", "CONTIG", "START", "STOP", "ORIENTATION", "CUT_OFF", "Best_Hit_evalue", "Best_Hit_ARO", "Best_Identites", "ARO", "ARO_name", "Model_type", "SNP", "AR0_category", "bit_score"])
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

			for it in data[item]:
				cgList = []
				if checkKeyExisted("ARO_category", data[item][it]):
					for aroctkey in data[item][it]["ARO_category"]:
						cgList.append(str(data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace')))

				if data[item][it]["model_type"] == "model-mutation":
					temp = data[item][it]["SNP"]["original"] + str(data[item][it]["SNP"]["position"]) + data[item][it]["SNP"]["change"]
					snpList.append(convert(temp))
				elif data[item][it]["model_type"] == "model-blastP":
					snpList.append("n/a")

				AROlist.append(convert(data[item][it]["ARO_accession"]))
				AROnameList.append(convert(data[item][it]["ARO_name"]))
				bitScoreList.append(data[item][it]["bit-score"])
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
				else:
					startCompare = True
					minevalue = data[item][it]["evalue"]
					minARO = data[item][it]["ARO_name"]

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
				writer.writerow([findnthbar(item, 0), findORFfrom(item), int(findnthbar(item, 4))-1, int(findnthbar(item, 5))-1, findnthbar(item, 3), ', '.join(list(clist)), minevalue, minARO, max(identityList), ', '.join(map(lambda x:"ARO:"+x, AROlist)), ', '.join(list(arocatset)), ', '.join(list(tl)), snpList, ', '.join(AROsortedList), ', '.join(map(str, bitScoreList))])


def main(afile):
	printCSV(afile)


if __name__ == '__main__':
	main(sys.argv[1])
	"""required: sys.argv[1] must be a json file"""
