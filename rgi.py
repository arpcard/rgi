import json
import sys
import os
import fqToFsa
import contigToProteins
#import blastnsnp


#remove temporary file
def removeTemp():
	'''
	if os.path.isfile("contigToPro.fasta"):
		os.remove("contigToPro.fasta")
	if os.path.isfile("read.fsa"):
		os.remove("read.fsa")
	if os.path.isfile("contig.fsa"):
		os.remove("contig.fsa")
	if os.path.isfile("contigToORF.fsa"):
		os.remove("contigToORF.fsa")
	if os.path.isfile("blastRes.xml"):
		os.remove("blastRes.xml")
	if os.path.isfile("blastpjson"):
		os.remove("blastpjson")
	if os.path.isfile("proteindb.fsa"):
		os.remove("proteindb.fsa")
	if os.path.isfile("dnadb.fsa"):
		os.remove("dnadb.fsa")
	if os.path.isfile("protein.db.phr"):
		os.remove("protein.db.phr")
		os.remove("protein.db.pin")
		os.remove("protein.db.psq")
	if os.path.isfile("temp.fsa"):
		os.remove("temp.fsa")
	if os.path.isfile("dna.db.nhr"):
		os.remove("dna.db.nhr")
		os.remove("dna.db.nin")
		os.remove("dna.db.nsq")
	'''

# check if a string is a float number
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


#output the information particular field from alignment.Title by splicing it by '|'
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


def checkBeforeBlast(inType, inputSeq):
	# give warning if inType is neither 'protein' nor 'dna' and then exit this program
	if inType != 'protein' and inType != 'contig' and inType != 'read':
		print>>sys.stderr, "inType must be one of protein, contig, orf or read. "
		removeTemp()
		exit()


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


def writeFASTAfromJson():
	noSeqList = []

	if os.path.isfile("proteindb.fsa") == False:
		with open("card.json") as json_file:
			json_data = json.load(json_file)
			with open ('proteindb.fsa', 'w') as wp:

				for item in json_data:

					if item.isdigit(): #get rid of __comment __timestamp etc
						# model_type: blastP only (pass_evalue)
						if json_data[item]["model_type_id"] == "40292":
							pass_eval = 1e-30
							if checkKeyExisted("model_param", json_data[item]):
								if checkKeyExisted("blastp_evalue", json_data[item]["model_param"]):
									pass_eval = json_data[item]["model_param"]["blastp_evalue"]["param_value"]
							
							if checkKeyExisted("model_sequences", json_data[item]):
								for seqkey in json_data[item]["model_sequences"]["sequence"]:
									print>>wp, ('>' + item + '_' + seqkey + " | model_type_id: 40292" + " | pass_evalue: " + str(pass_eval))
									print>>wp, (json_data[item]["model_sequences"]["sequence"][seqkey]["protein_sequence"]["sequence"])

							else:
								 noSeqList.append(item)

						# model_type: blastP + SNP (pass_evalue + snp)
						elif json_data[item]["model_type_id"] == "40293":
							snpList = ""
							if checkKeyExisted('snp', json_data[item]['model_param']):
								for key in json_data[item]['model_param']['snp']['param_value']:
									snpList += json_data[item]['model_param']['snp']['param_value'][key]
									snpList += ','
							pass_eval = 1e-30
							if checkKeyExisted("model_param", json_data[item]):
								if checkKeyExisted("blastp_evalue", json_data[item]["model_param"]):
									pass_eval = json_data[item]["model_param"]["blastp_evalue"]["param_value"]
							
							if checkKeyExisted("model_sequences", json_data[item]):
								for seqkey in json_data[item]["model_sequences"]["sequence"]:
									print>>wp, ('>' + item + '_' + seqkey + " | model_type_id: 40293" + " | pass_evalue: " + str(pass_eval) + " | SNP: " + snpList)
									print>>wp, (json_data[item]["model_sequences"]["sequence"][seqkey]["protein_sequence"]["sequence"])
							else:
								 noSeqList.append(item)
		json_file.close()
	#get a list of models who are incomplete
	'''if noSeqList:
		print>>sys.stderr, noSeqList'''		


#make protein (and dna if needed) database:
def makeBlastDB(inType, inputSeq):
	if os.path.isfile("proteindb.fsa") == True and os.path.exists("proteindb.fsa") == True  and os.path.exists("protein.db.phr") == True and os.path.exists("protein.db.pin") == True and os.path.exists("protein.db.psq") == True :
		print "DB exists"
	else:
		print "create DB."
		os.system('makeblastdb -in proteindb.fsa -dbtype prot -out protein.db')
	#os.system('makeblastdb -in proteindb.fsa -dbtype prot -out protein.db')

def runBlast(inType, inputSeq):	
	startBlast = False
	if inType == 'contig':
		contigToProteins.main(inputSeq)
		if os.stat("contig.fsa").st_size != 0:
			from Bio.Blast.Applications import NcbiblastpCommandline
			blastCLine = NcbiblastpCommandline(query="./contig.fsa", db="protein.db", outfmt=5, out="blastRes.xml")
			stdt, stdr = blastCLine()
			result_handle = open("blastRes.xml")
			startBlast = True
		
	elif inType == 'protein':
		from Bio.Blast.Applications import NcbiblastpCommandline
		blastCLine = NcbiblastpCommandline(query=inputSeq, db="protein.db", outfmt=5, out="blastRes.xml")
		stdt, stdr = blastCLine()
		result_handle = open("blastRes.xml")
		startBlast = True

	elif inType == 'read':
		fqToFsa.main(inputSeq)
		from Bio.Blast.Applications import NcbiblastxCommandline
		blastCLine = NcbiblastxCommandline(query='./read.fsa', db="protein.db", outfmt=5, out="blastRes.xml")
		stdt, stdr = blastCLine()
		result_handle = open("blastRes.xml")
		startBlast = True
	
	if startBlast:
		from Bio.Blast import NCBIXML
		blast_records = NCBIXML.parse(result_handle)
		blastResults = {}
		pjson = {}


		with open("card.json") as json_file:
			json_data = json.load(json_file)

		for blast_record in blast_records:
			perfect = {}
			strict = {}
			loose = {}


			for alignment in blast_record.alignments:			
				alignTitle = alignment.title
				orfInfo = blast_record.query.encode('ascii','replace')
				c = 0
				barc = 0
				for eachc in orfInfo:
					if barc >= 6:
						break
					elif eachc == '|':
						barc += 1
						c += 1
					else:
						c += 1
				orffrom = orfInfo[c:]

				modelTypeID = findnthbar(alignTitle, 0)
				spacepos = alignTitle.index(' ')
				hitid = alignTitle[0:spacepos]
				hitid = hitid.encode('ascii','replace')

				modelDescrpt =alignTitle[alignTitle.index(' ')+1:]
				underscoreinMD = modelDescrpt.index('_')
				modelID = modelDescrpt[0:underscoreinMD]
				seqinModel = modelDescrpt[underscoreinMD+1: modelDescrpt.index(' ')]

				modelTypeDscp = alignTitle[alignTitle.index(':')+2:]
				modelTypeId = modelTypeDscp[0:modelTypeDscp.index(' ')]
				passevalue = modelTypeDscp [modelTypeDscp.index(':')+2:]

				if isfloat(passevalue):
					truePassEvalue = passevalue
				else:
					truePassEvalue = float(passevalue[0:passevalue.find(' ')])

				if modelTypeID == 40293:
					init = 0
					evalueSNP = findnthbar(alignTitle, 2)
					snpL = []
					snpdictlist = []
					temp = ""

					for eachc in evalueSNP:
						if eachc == ',':
							snpL.append(temp)
							temp = ""
						else:
							temp += eachc

					for eachsnp in snpL:					
						snpdictlist.append({"original": eachsnp[0], "change": eachsnp[-1], "position": int(eachsnp[1:-1])})

					for hsp in alignment.hsps:
						'''print>>sys.stderr, hsp.identities
						print>>sys.stderr, len(hsp.match)
						print "hsp.query: \n" +hsp.query + "\n"
						print "hsp.sbj: \n" +hsp.sbjct + "\n"'''

						querySeq =  hsp.query.replace('-', '')
						realQueryLength = len(querySeq) 

						sbjctSeq = hsp.sbjct.replace('-', '') 
						realSbjctLength = len(sbjctSeq) 
						sinsidedict = {}

						'''print "QuerySeq: \n" +querySeq + "\n"
						print "sbjctSeq: \n" +sbjctSeq + "\n"
						print "realQueryLength: \n" + str(realQueryLength) + "\n"
						print "realSbjctLength: \n" +str(realSbjctLength) + "\n"
						print "subject_start: \n" + str(hsp.sbjct_start) + "\n"
						print "query_start: \n" + str(hsp.query_start) + "\n"'''
						


							
						for eachs in snpdictlist:
							pos = eachs["position"]
							ori = eachs["original"]
							chan = eachs["change"]

							if hsp.sbjct_start < pos and (hsp.sbjct_start + realSbjctLength) > pos:
								'''print>>sys.stderr, eachs
								print>>sys.stderr, hsp.sbjct
								print>>sys.stderr, hsp.sbjct_start
								print>>sys.stderr, ("wrong snp in query: " + hsp.query[pos - hsp.sbjct_start])
								print>>sys.stderr, hsp.query'''
								'''c = 0
								target = pos - hsp.sbjct_start
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
										c += 1'''
								'''print>>sys.stderr, ("corret: snp in sbject: " + sbjctSeq[target])
								print>>sys.stderr, ("correct: snp in query: " + hsp.query[snpInQuery])'''

								#if hsp.query[snpInQuery] == chan and sbjctSeq[target] == ori:
								if hsp.query[pos - hsp.sbjct_start + findNumDash(hsp.sbjct, (pos-hsp.sbjct_start))] == chan and hsp.sbjct[pos - hsp.sbjct_start +findNumDash(hsp.sbjct, (pos-hsp.sbjct_start))] == ori:								# 224 = pos. sbject_start = 9 = 216==

								#pos = 224, start = 9 = 216
									if hsp.expect <= truePassEvalue:
										sinsidedict = {}
										sinsidedict["type_match"] = "Strict"
										sinsidedict["SNP"] = eachs
										sinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
										sinsidedict["orf_start"] = findnthbar(orfInfo, 1)
										sinsidedict["orf_end"] = findnthbar(orfInfo, 2)
										sinsidedict["orf_From"] = orffrom
										
										sinsidedict["model_name"] = json_data[modelID]["model_name"]
										sinsidedict["model_type"] = json_data[modelID]["model_type"]
										sinsidedict["model_id"] = modelID
										sinsidedict["pass_evalue"] = passevalue
										sinsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
										sinsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
										sinsidedict["ARO_category"] = json_data[modelID]["ARO_category"]
										sinsidedict["evalue"] = hsp.expect
										sinsidedict["max-identities"] = hsp.identities
										sinsidedict["bit-score"] = hsp.bits
										if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
											sinsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]

										sinsidedict["query"] = hsp.query.encode('ascii','replace')
										sinsidedict["match"] = hsp.match.encode('ascii','replace')
										sinsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
										sinsidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]

										if inType == 'contig':
											sinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
											sinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
											sinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
											sinsidedict["orf_start"] = findnthbar(orfInfo, 1)
											sinsidedict["orf_end"] = findnthbar(orfInfo, 2)
											sinsidedict["orf_From"] = orffrom

										elif inType == 'protein':
											sinsidedict["query_start"] = hsp.query_start
											sinsidedict["query_end"] = hsp.query_start + realQueryLength
											sinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
												
										strict[hitid + "|hsp_num:" + str(init)] = sinsidedict
										init += 1

									else:
										slinsidedict = {}
										slinsidedict["type_match"] = "Loose"
										slinsidedict["SNP"] = eachs
										slinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
										slinsidedict["orf_start"] = findnthbar(orfInfo, 1)					
										slinsidedict["orf_end"] = findnthbar(orfInfo, 2)
										slinsidedict["orf_From"] = orffrom

										slinsidedict["model_name"] = json_data[modelID]["model_name"]
										slinsidedict["model_type"] = json_data[modelID]["model_type"]
										slinsidedict["pass_evalue"] = passevalue
										slinsidedict["model_id"] = modelID
										slinsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
										slinsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
										slinsidedict["evalue"] = hsp.expect
										slinsidedict["bit-score"] = hsp.bits
										slinsidedict["max-identities"] = hsp.identities
										if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
											slinsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]

										slinsidedict["query"] = hsp.query.encode('ascii','replace')
										slinsidedict["match"] = hsp.match.encode('ascii','replace')
										slinsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
										slinsidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]

										if inType == 'contig':
											slinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
											slinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
											slinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
											slinsidedict["orf_start"] = findnthbar(orfInfo, 1)
											slinsidedict["orf_end"] = findnthbar(orfInfo, 2)
											slinsidedict["orf_From"] = orffrom

										elif inType == 'protein':
											slinsidedict["query_start"] = hsp.query_start
											slinsidedict["query_end"] = hsp.query_start + realQueryLength
											slinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
												
										loose[hitid + "|hsp_num:" + str(init)] = slinsidedict
										init += 1
		

				elif modelTypeID == 40292:
					init = 0
					passevalue = findnthbar(alignTitle, 1)

					for hsp in alignment.hsps:
						'''print>>sys.stderr, hsp.identities
						print>>sys.stderr, alignment.length
						print>>sys.stderr, hsp.query
						print>>sys.stderr, hsp.query[hsp.query_start]'''

						querySeq = hsp.query.replace('-', '')
						realQueryLength = len(querySeq)

						if hsp.query == json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]:
							ppinsidedict = {}
							ppinsidedict["type_match"] = "Perfect"
							ppinsidedict["model_id"] = modelID
							ppinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
							ppinsidedict["orf_start"] = findnthbar(orfInfo, 1)
							ppinsidedict["orf_end"] = findnthbar(orfInfo, 2)
							ppinsidedict["orf_From"] = orffrom

							ppinsidedict["model_name"] = json_data[modelID]["model_name"]
							ppinsidedict["model_type"] = json_data[modelID]["model_type"]
							ppinsidedict["pass_evalue"] = passevalue
							ppinsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
							ppinsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
							ppinsidedict["ARO_category"] = json_data[modelID]["ARO_category"]
							ppinsidedict["evalue"] = hsp.expect
							ppinsidedict["bit-score"] = hsp.bits
							ppinsidedict["max-identities"] = hsp.identities
							if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
								ppinsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
							
							ppinsidedict["query"] = hsp.query.encode('ascii','replace')
							ppinsidedict["match"] = hsp.match.encode('ascii','replace')
							ppinsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
							ppinsidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]

							if inType == 'contig':
								ppinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
								ppinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
								ppinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
								ppinsidedict["orf_start"] = findnthbar(orfInfo, 1)
								ppinsidedict["orf_end"] = findnthbar(orfInfo, 2)
								ppinsidedict["orf_From"] = orffrom

							elif inType == 'protein':
								ppinsidedict["query_start"] = hsp.query_start
								ppinsidedict["query_end"] = hsp.query_start + realQueryLength
								ppinsidedict["query_From"] = blast_record.query.encode('ascii','replace')

							perfect[hitid + "|hsp_num:" + str(init)] = ppinsidedict
							init += 1
									
						elif hsp.expect <= passevalue:
							insidedict = {}
							insidedict["type_match"] = "Strict"
							insidedict["orf_strand"] = findnthbar(orfInfo, 0)
							insidedict["orf_start"] = findnthbar(orfInfo, 1)							
							insidedict["orf_end"] = findnthbar(orfInfo, 2)
							insidedict["orf_From"] = orffrom

							insidedict["model_name"] = json_data[modelID]["model_name"]
							insidedict["model_type"] = json_data[modelID]["model_type"]
							insidedict["model_id"] = modelID
							insidedict["pass_evalue"] = passevalue
							insidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
							insidedict["ARO_name"] = json_data[modelID]["ARO_name"]
							insidedict["ARO_category"] = json_data[modelID]["ARO_category"]
							insidedict["evalue"] = hsp.expect
							insidedict["bit-score"] = hsp.bits
							insidedict["max-identities"] = hsp.identities
							if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
								insidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
							
							insidedict["query"] = hsp.query.encode('ascii','replace')
							insidedict["match"] = hsp.match.encode('ascii','replace')
							insidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
							insidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]

							if inType == 'contig':
								insidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
								insidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
								insidedict["orf_strand"] = findnthbar(orfInfo, 0)
								insidedict["orf_start"] = findnthbar(orfInfo, 1)
								insidedict["orf_end"] = findnthbar(orfInfo, 2)
								insidedict["orf_From"] = orffrom

							elif inType == 'protein':
								insidedict["query_start"] = hsp.query_start
								insidedict["query_end"] = hsp.query_start + realQueryLength
								insidedict["query_From"] = blast_record.query.encode('ascii','replace')

							strict[hitid + "|hsp_num:" + str(init)] = insidedict
							init += 1
							
						else:
							linsidedict = {}
							linsidedict["type_match"] = "Loose"
							linsidedict["orf_strand"] = findnthbar(orfInfo, 0)
							linsidedict["orf_start"] = findnthbar(orfInfo, 1)
							linsidedict["orf_end"] = findnthbar(orfInfo, 2)
							linsidedict["orf_From"] = orffrom
							
							linsidedict["model_name"] = json_data[modelID]["model_name"]
							linsidedict["model_type"] = json_data[modelID]["model_type"]
							linsidedict["pass_evalue"] = passevalue
							linsidedict["model_id"] = modelID
							linsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
							linsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
							linsidedict["ARO_category"] = json_data[modelID]["ARO_category"]
							linsidedict["evalue"] = hsp.expect
							linsidedict["max-identities"] = hsp.identities
							linsidedict["bit-score"] = hsp.bits
							if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
								linsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]

							linsidedict["query"] = hsp.query.encode('ascii','replace')
							linsidedict["match"] = hsp.match.encode('ascii','replace')
							linsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
							linsidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]

							if inType == 'contig':
								linsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
								linsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
								linsidedict["orf_strand"] = findnthbar(orfInfo, 0)
								linsidedict["orf_start"] = findnthbar(orfInfo, 1)
								linsidedict["orf_end"] = findnthbar(orfInfo, 2)
								linsidedict["orf_From"] = orffrom

							elif inType == 'protein':
								linsidedict["query_start"] = hsp.query_start
								linsidedict["query_end"] = hsp.query_start + realQueryLength
								linsidedict["query_From"] = blast_record.query.encode('ascii','replace')
									
							loose[hitid + "|hsp_num:" + str(init)] = linsidedict
							init += 1
									
			if len(perfect) == 0 and len(strict) == 0:
				blastResults[blast_record.query.encode('ascii','replace')] = loose
				pscore = 1
				
			elif len(perfect) == 0:
				blastResults[blast_record.query.encode('ascii','replace')] = strict
				pscore = 2
				
			else:
				blastResults[blast_record.query.encode('ascii','replace')] = perfect
				pscore = 3

		pjson = json.dumps(blastResults)

		with open("Report.json", 'w') as f:
			print>>f, pjson
		if inType == 'contig':
			for gene in blastResults:
				for hsp in blastResults[gene]:
					if blastResults[gene][hsp]["orf_end"] < blastResults[gene][hsp]["query_end"]:
						print>>sys.stderr, hsp
						print>>sys.stderr, blastResults[gene][hsp]
	return pjson


def findNumDash(subject, index):
	numDash = 0
	toFind = "-"
	stringList = list(subject)
	output = []
	
	for i in range(0, len(stringList)):
		if (stringList[i] == toFind): 
			numDash += 1
		else: 
			output.append(stringList[i])
		if (len(output) == index): 
			break
	return numDash
	


"""
function get_number_dash($subject,$index){

    $number_of_dash = 0;
    $toFind = "-";
    $array = str_split($subject);

    $output = [];

    for($i = 0; $i < strlen($subject); $i++){

        if($array[$i] == $toFind){
            $number_of_dash++;
        }else{
            $output[] = $array[$i];
        }

        if(count($output) == $index){ break; }
    }

    return $number_of_dash;
}
"""

def main(inType, inputSeq):

	checkBeforeBlast(inType, inputSeq)
	writeFASTAfromJson()
	makeBlastDB(inType, inputSeq)

	try:	
		if inType == 'contig':
			bpjson = runBlast(inType, inputSeq)
			'''#1. uncomment here to add blastn+snp model to RGI if we have extended our database, 
			#2. uncomment the import line above,
			#3. return myjson instead of bpjson in this case
			njson = blastnsnp.main('contig', inputSeq)
			myjson = dict(bpjson.item() + njson.item())
			removeTemp()
			return myjson'''
			removeTemp()
			return bpjson
		else:
			bpjson = runBlast(inType, inputSeq)
			removeTemp()
			return bpjson

	except Exception as inst:
		#pass
		print>>sys.stderr, inst
		removeTemp()



if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
	"""required: 
		argv[1] must be one of contig, orf, protein, read
		inputSeq must be in either FASTA (contig and protein) or FASTQ(read) format!"""
