import sys
import json
import os


def removeTempFile():	
	if os.path.isfile("blastnRes.xml"):
		os.remove("blastnRes.xml")
	if os.path.isfile("dnadb.fsa"):
		os.remove("dnadb.fsa")
	if os.path.isfile("dna.db.nhr"):
		os.remove("dna.db.nhr")
		os.remove("dna.db.nin")
		os.remove("dna.db.nsq")


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


def convert(input):
	if isinstance(input, dict):
		return dict((convert(key), convert(value)) for key, value in input.iteritems())
	elif isinstance(input, list):
		return [convert(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input


def findnthbar(bunchstr, n):
	barc = 0
	start = n+3
	over = n+4
	temp = ""
	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if barc == start:
			temp += eachc
		if barc == over:
			break
	colonpos = temp.find(':')
	res=temp[colonpos+2:]
	res=res.rstrip()
	if res.isdigit():
		return int(res)
	elif isfloat(res):
		return float(res)
	else:
		res = res.encode('ascii','replace')
		return res


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


def writeFASTAfromJson():
	if os.path.isfile("dnadb.fsa") == False:
		with open("card.json") as json_file:
			json_data = json.load(json_file)
			
			with open ('dnadb.fsa', 'w') as wd:
				for item in json_data:
					if item.isdigit():
						# model_type: blastN + SNP (pass_evalue + snp)
						if json_data[item]["model_type_id"] == "40295":
							snpList = ""
							if checkKeyExisted('snp', json_data[item]['model_param']):
								for key in json_data[item]['model_param']['snp']['param_value']:
									snpList += json_data[item]['model_param']['snp']['param_value'][key]
									snpList += ','

							if checkKeyExisted("model_sequences", json_data[item]):
								for seqkey in json_data[item]["model_sequences"]["sequence"]:
									print>>wd, ('>' + item + '_' + seqkey + " | model_type_id: 40295" + " | SNP: " + snpList)
									print>>wd, (json_data[item]["model_sequences"]["sequence"][seqkey]["dna_sequence"]["sequence"])
			wd.close()
		json_file.close()
		

def runBlastnSnp (inType, inputSeq):
	if os.path.isfile("dna.db.nsq") == False: 
		os.system('makeblastdb -in dnadb.fsa -dbtype nucl -out dna.db')		
	os.system("blastn -out blastnRes.xml -outfmt 5 -query " + inputSeq + " -db dna.db")
	
	result_handle = open('blastnRes.xml')
	from Bio.Blast import NCBIXML
	blast_records = NCBIXML.parse(result_handle)			
	blastnResults = {}
	result_handle.close()

	with open("card.json") as json_file:
		json_data = json.load(json_file)
	json_file.close()

	for blast_record in blast_records:
		nsnp = {}

		for alignment in blast_record.alignments:
			alignTitle = alignment.title
			modelTypeID = findnthbar(alignTitle, 0)
			spacepos = alignTitle.index(' ')
			hitid = alignTitle[0:spacepos]
			hitid = hitid.encode('ascii','replace')

			modelDescrpt = alignTitle[alignTitle.index(' ')+1:]
			underscoreinMD = modelDescrpt.index('_')
			modelID = modelDescrpt[0:underscoreinMD]
			seqinModel = modelDescrpt[underscoreinMD+1:modelDescrpt.index(' ')]

			modelTypeDscp = alignTitle[alignTitle.index(':')+2:]
			modelTypeId = modelTypeDscp[0:modelTypeDscp.index(' ')]
			passevalue = modelTypeDscp [modelTypeDscp.index(':')+2:]

			init = 0
			evalueSNP = findnthbar(alignTitle, 1)
				#print evalueSNP
			snpL = []
			temp = ""

			for eachc in evalueSNP:
				if eachc == ',':
					snpL.append(temp)
					temp = ""
				else:
					temp += eachc

			snpdictlist = []
			for eachsnp in snpL:					
				snpdictlist.append({"original": eachsnp[0], "change": eachsnp[-1], "position": int(eachsnp[1:-1])})
					
			for hsp in alignment.hsps:
				querySeq = hsp.query.replace('-', '')
				realQueryLength = len(querySeq)

				sbjctSeq = hsp.sbjct.replace('-', '')
				realSbjctLength = len(sbjctSeq)
				sinsidedict = {}

				for eachs in snpdictlist:
					pos = eachs["position"]
					ori = eachs["original"]
					chan = eachs["change"]

					if hsp.sbjct_start <= pos and (hsp.sbjct_start + realSbjctLength) > pos:
						target = pos - hsp.sbjct_start
						c = 0
						snpInQuery = 0

						for ch in hsp.sbjct:
							if c == target:
								if ch == '-':
									snpInQuery += 1
								else:
									break
							elif ch == '-':
								snpInQuery += 1
							else:
								snpInQuery += 1
								c += 1

						if hsp.query[snpInQuery] == chan and sbjctSeq[target] == ori:
							nsinsidedict = {}
							nsinsidedict["model_id"] = modelID
							nsinsidedict["SNP"] = eachs
							nsinsidedict["query_start"] = hsp.query_start
							nsinsidedict["query_end"] = hsp.query_start + realQueryLength - 1
							nsinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
							
							if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
								nsinsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
							nsinsidedict["model_name"] = json_data[modelID]["model_name"]
							nsinsidedict["model_type"] = json_data[modelID]["model_type"]
							nsinsidedict["pass_evalue"] = passevalue
							nsinsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
							nsinsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
							nsinsidedict["ARO_category"] = json_data[modelID]["ARO_category"]
							nsinsidedict["evalue"] = hsp.expect
							nsinsidedict["max-identites"] = hsp.identities
							nsinsidedict["bit-score"] = hsp.bits
							
							nsinsidedict["query"] = hsp.query.encode('ascii','replace')
							nsinsidedict["match"] = hsp.match.encode('ascii','replace')
							nsinsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
												
							nsnp[hitid + "|hsp_num:" + str(init)] = nsinsidedict
							init += 1

			blastnResults[blast_record.query.encode('ascii','replace')] = nsnp	
	jsonRes = json.dumps(blastnResults)
	with open ('blastnReport.json', 'w') as wj:
		print>>wj, jsonRes
	wj.close()
	return jsonRes


def main(inType, inputSeq):
	writeFASTAfromJson()
	try:
		myjson = runBlastnSnp('contig', inputSeq)
		removeTempFile()
		return myjson

	except Exception as inst:
		#pass
		print>>sys.stderr, inst
		removeTempFile()


if __name__ == '__main__':
	main('contig', sys.argv[2])
	"""required: inputSeq must be in FASTA format!"""
