import json
import sys
import os
import fqToFsa
import contigToProteins
import contigToORF
#import blastnsnp
import argparse
import filepaths
import gzip, zlib
#from gzopen import gzopen
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import logging
import hashlib
import time
import convertJsonToTSV

script_path = filepaths.determine_path()
working_directory = os.getcwd()

clean_files = []

path = script_path+"/_db/"
data_path = script_path+"/_data/"
tmp = script_path+"/_tmp/"

#remove temporary file
def removeTemp():
	for tempfile in clean_files:
		if os.path.isfile(tempfile):
			logging.info("removeTemp => remove temp file: " + str(tempfile))
			os.remove(tempfile)


# check if a string is a float number
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


#output the information particular field from alignment.Title by splicing it by '|'
def _findnthbar(bunchstr, n):
	if bunchstr[0:5] in "gene_":
		#print bunchstr[0:5]
		pass
	else:
		bunchstr = "_|_|_ "+ bunchstr

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

def findnthbar(bunchstr, n):
	if bunchstr[0:5] in "gene_":
		pass
	else:
		if "model_type_id" in bunchstr:
			pass
		else:
			bunchstr = "_|_|_ "+ bunchstr

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

#output the information particular field from alignment.Title by splicing it by '#'
def findnthbar2(bunchstr, n):

	if "#" in bunchstr:
		arr = bunchstr.split("#")
		if n < len(arr):
			# gene id
			if n == 1:
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
	else:
		return 0

def checkBeforeBlast(inType, inputSeq):
	logging.info("checkBeforeBlast => " + str([inType,inputSeq]))
	# give warning if inType is neither 'protein' nor 'dna' and then exit this program
	if inType != 'protein' and inType != 'contig' and inType != 'read':
		print>>sys.stderr, "inType must be one of protein, contig, orf or read. "
		logging.error("checkBeforeBlast => " + str(inType) + ", inType must be one of protein, contig, orf or read. ")
		removeTemp()
		exit()
	logging.info("checkBeforeBlast => success")


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


def writeFASTAfromJson():
	noSeqList = []

	if os.path.isfile(path+"proteindb.fsa") == False:
		with open(data_path+"card.json") as json_file:
			json_data = json.load(json_file)
			with open (path+'proteindb.fsa', 'w') as wp:

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
							pass_eval = 1e-30
							if checkKeyExisted("model_param", json_data[item]):
								if checkKeyExisted("blastp_evalue", json_data[item]["model_param"]):
									pass_eval = json_data[item]["model_param"]["blastp_evalue"]["param_value"]
									
								if checkKeyExisted("snp", json_data[item]['model_param']):
									for key in json_data[item]['model_param']['snp']['param_value']:
										snpList += json_data[item]['model_param']['snp']['param_value'][key]
										snpList += ','
							
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
	if os.path.isfile(path+"proteindb.fsa") == True and os.path.exists(path+"proteindb.fsa") == True  and os.path.exists(path+"protein.db.phr") == True and os.path.exists(path+"protein.db.pin") == True and os.path.exists(path+"protein.db.psq") == True :
		print "blast DB exists"
	else:
		print "create blast DB."
		os.system('makeblastdb -in '+path+'proteindb.fsa -dbtype prot -out '+path+'protein.db')
	#os.system('makeblastdb -in proteindb.fsa -dbtype prot -out protein.db')

def makeDiamondDB():
	if os.path.isfile(path+"proteindb.fsa") == True and os.path.exists(path+"proteindb.fsa") == True  and os.path.exists(path+"protein.db.dmnd") == True:
		print "diamond DB exists"
	else:
		print "create diamond DB."
		os.system('diamond makedb --in '+path+'proteindb.fsa --db '+path+'protein.db')

def getORFDNASequence(file_name,orf):
	predicted_genes_dict = {}
	if os.stat(working_directory+"/"+file_name+".contigToORF.fsa").st_size != 0:
		from Bio import SeqIO
		clean_files.append(working_directory+"/"+file_name+".contigToORF.fsa")
		for record in SeqIO.parse(working_directory+"/"+file_name+".contigToORF.fsa", 'fasta'):
			if orf == "genemark" and '|' in record.id:
				predicted_genes_dict[record.id[:str(record.id).index('|')]] = str(record.seq)
			else:
				predicted_genes_dict[record.id] = str(record.seq)
	
	# write json for all predicted file
	pjson = json.dumps(predicted_genes_dict)
	clean_files.append(working_directory+'/'+file_name+'.predictedGenes.json')
	with open(working_directory+'/'+file_name+'.predictedGenes.json', 'w') as wf:
		print>>wf, pjson
	wf.close()

	return predicted_genes_dict

def getSubmittedProteinSequence(afile):
	submitted_proteins_dict = {}
	if os.stat(afile).st_size != 0:
		from Bio import SeqIO
		for record in SeqIO.parse(afile, 'fasta'):
			submitted_proteins_dict[record.id] = str(record.seq)

	return submitted_proteins_dict

def runDiamond(inType, inputSeq, threads, outputFile):

	cfilter = "	--index-chunks 1 --block-size 1 --evalue 10" #--quiet

	if inType == 'contig' or inType == 'read':
		logging.info("runDiamond => blastx")
		logging.info('runDiamond => diamond blastx --in '+path+'proteindb.fsa --db '+path+'protein.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles --sensitive'+cfilter)
		os.system('diamond blastx --in '+path+'proteindb.fsa --db '+path+'protein.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles --more-sensitive'+cfilter)
	else:
		logging.info("runDiamond => blastp")
		logging.info('runDiamond => diamond blastp --in '+path+'proteindb.fsa --db '+path+'protein.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles'+cfilter)
		os.system('diamond blastp --in '+path+'proteindb.fsa --db '+path+'protein.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles'+cfilter)
    
	logging.info('runDiamond => diamond view --outfmt xml --daa '+inputSeq+'.daa --out '+outputFile + cfilter)
	os.system('diamond view --outfmt xml --daa '+inputSeq+'.daa --out '+outputFile + cfilter )
	logging.info('runDiamond => diamond view --daa '+inputSeq+'.daa --out '+outputFile+".txt" + cfilter)
	os.system('diamond view --daa '+inputSeq+'.daa --out '+outputFile+".txt" + cfilter )

	clean_files.append(inputSeq+'.daa')
	clean_files.append(outputFile+'.txt')

def getHashName(name):
	m = hashlib.md5()
	t = time.gmtime()
	m.update(name + str(t))
	return m.hexdigest()

def encodeHeaders(fastaFile):
	_fastaFile = os.path.basename(fastaFile)
	headers_dict = {}
	ofile = open("my_fasta.txt", "w")
	from Bio import SeqIO
	for record in SeqIO.parse(fastaFile, 'fasta'):
		new_header = getHashName(record.id)
		headers_dict[new_header] = record.id
		ofile.write(">" + new_header + "\n" + str(record.seq) + "\n")
	ofile.close()
	os.rename("my_fasta.txt", fastaFile)
	# write dictionaty for encoded headers
	headers_dict_json = json.dumps(headers_dict)
	clean_files.append(working_directory+"/"+_fastaFile+".encodedHeaders.json")
	with open(working_directory+"/"+_fastaFile+".encodedHeaders.json", 'w') as f:
	 	print>>f, headers_dict_json
	f.close()	



def runBlast(inType, inputSeq, threads, outputFile, criteria, data_type, clean, orf, alignment_tool):	

	pjson = {}
	startBlast = False
	predicted_genes_dict = {}
	submitted_proteins_dict = {}

	# encode header using md5 and store this to dictionary - this will solve issues with headers which have same characters as alignment tools :)
	#encodeHeaders(inputSeq)
	file_name = os.path.basename(inputSeq)
	clean_files.append(working_directory+"/"+file_name+".blastRes.xml")
	#print ">>>>", working_directory+"/"+file_name+".blastRes.xml"

	if inType == 'contig':
		if orf == "genemark":
			logging.info("runBlast => contigToProteins => start")
			contigToProteins.main(inputSeq,clean)
			logging.info("runBlast => contigToProteins => done")

			# get predicted dna
			logging.info("runBlast => contigToORF => start")
			contigToORF.main(inputSeq,clean,orf)
			logging.info("runBlast => contigToORF => done")

			logging.info("runBlast => getORFDNASequence => start")
			predicted_genes_dict = getORFDNASequence(file_name,orf)
			logging.info("runBlast => getORFDNASequence => done")
		else:
			#logging.info("runBlast => contigToProteins => start")
			#contigToProteins.main(inputSeq,clean)
			#logging.info("runBlast => contigToProteins => done")

			# get predicted dna
			logging.info("runBlast => contigToORF => start")
			contigToORF.main(inputSeq,clean,orf)
			logging.info("runBlast => contigToORF => done")

			logging.info("runBlast => getORFDNASequence => start")
			predicted_genes_dict = getORFDNASequence(file_name,orf)
			logging.info("runBlast => getORFDNASequence => done")			

		if os.stat(working_directory+"/"+file_name+".contig.fsa").st_size != 0:

			logging.info("runBlast => start blastP for inType: " + inType)
			clean_files.append(working_directory+"/"+file_name+".contig.fsa")

			if alignment_tool == "diamond":
				runDiamond("protein", working_directory+"/"+file_name+".contig.fsa", threads, working_directory+"/"+file_name+".blastRes.xml")
			else:
				from Bio.Blast.Applications import NcbiblastpCommandline
				blastCLine = NcbiblastpCommandline(query=working_directory+"/"+file_name+".contig.fsa", db=path+"protein.db", outfmt=5, out=working_directory+"/"+file_name+".blastRes.xml",num_threads=threads)#,num_alignments=1)
				stdt, stdr = blastCLine()

			result_handle = open(working_directory+"/"+file_name+".blastRes.xml")
			startBlast = True

	elif inType == 'protein':

		submitted_proteins_dict = getSubmittedProteinSequence(inputSeq)
		logging.info("runBlast => start blastP for inType: " + inType)

		if alignment_tool == "diamond":
			runDiamond(inType, inputSeq, threads, working_directory+"/"+file_name+".blastRes.xml")
		else:
			from Bio.Blast.Applications import NcbiblastpCommandline
			blastCLine = NcbiblastpCommandline(query=inputSeq, db=path+"protein.db", outfmt=5, out=working_directory+"/"+file_name+".blastRes.xml",num_threads=threads,num_alignments=1)
			stdt, stdr = blastCLine()

		result_handle = open(working_directory+"/"+file_name+".blastRes.xml")
		startBlast = True

	elif inType == 'read':
		logging.info("runBlast => fqToFsa => start")
		fqToFsa.main(inputSeq)
		logging.info("runBlast => fqToFsa => done")
		clean_files.append(working_directory+"/"+file_name+".read.fsa")
		logging.info("runBlast => start blastX for inType: " + inType)

		if alignment_tool == "diamond":
			runDiamond(inType, working_directory+"/"+file_name+".read.fsa", threads, working_directory+"/"+file_name+".blastRes.xml")
		else:
			from Bio.Blast.Applications import NcbiblastxCommandline
			blastCLine = NcbiblastxCommandline(query=working_directory+"/"+file_name+".read.fsa", db=path+"protein.db", outfmt=5, out=working_directory+"/"+file_name+".blastRes.xml",num_threads=threads,num_alignments=1)
			stdt, stdr = blastCLine()
		result_handle = open(working_directory+"/"+file_name+".blastRes.xml")
		startBlast = True

	if os.path.isfile(working_directory+"/"+file_name+".blastRes.xml"):
		# Check if we have any results
		if os.stat(working_directory+"/"+file_name+".blastRes.xml").st_size == 0:
			return pjson
	else:
		print "No results"
		return pjson

	if startBlast:
		logging.info("runBlast => startBlast = " + str(startBlast))
		from Bio.Blast import NCBIXML
		blast_records = NCBIXML.parse(result_handle)
		blastResults = {}
		pjson = {}

		with open(data_path+"card.json") as json_file:
			json_data = json.load(json_file)
		json_file.close()

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

				logging.info("runBlast => [info] | modelTypeID = " + str(alignTitle))

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

								card_sequence = str(json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"])
								orf_protein_sequence = ""
 
								if orf == "genemark":
									if predicted_genes_dict:
										orf_protein_sequence = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
									if submitted_proteins_dict:
										orf_protein_sequence = str(submitted_proteins_dict[orfInfo.split(" ")[0]])
								else:
									if predicted_genes_dict:
										orf_protein_sequence = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
									if submitted_proteins_dict:
										orf_protein_sequence = str(submitted_proteins_dict[orfInfo.split(" ")[0]])									

								logging.info("runBlast => [info] | Model:"+str(modelID) + " pos:" +str(pos) +" | "+str(hsp.query[pos - hsp.sbjct_start + findNumDash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND " + str(hsp.sbjct[pos - hsp.sbjct_start +findNumDash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

								# Report ONLY if the SNPs are present
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
										sinsidedict["model_type_id"] = modelTypeID

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
										sinsidedict["dnaSequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["dna_sequence"]["sequence"]
										
										if inType == 'contig':
											if orf == "genemark":
												sinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
												sinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
												sinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
												sinsidedict["orf_start"] = findnthbar(orfInfo, 1)
												sinsidedict["orf_end"] = findnthbar(orfInfo, 2)
												sinsidedict["orf_From"] = orffrom

												if orfInfo[:orfInfo.index('|')] in predicted_genes_dict:
													sinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('|')]] 
													sinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
												else:
													sinsidedict["orf_dna_sequence"] = ""
													sinsidedict["orf_prot_sequence"] = ""
											else:
												sinsidedict["query_start"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3
												sinsidedict["query_end"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
												sinsidedict["orf_strand"] = findnthbar2(orfInfo, 3)
												sinsidedict["orf_start"] = findnthbar2(orfInfo, 1)
												sinsidedict["orf_end"] = findnthbar2(orfInfo, 2)
												sinsidedict["orf_From"] = findnthbar2(orfInfo, 0)

												if orfInfo[:orfInfo.index('#')] in predicted_genes_dict:
													sinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('#')]] 
													sinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
												else:
													sinsidedict["orf_dna_sequence"] = ""
													sinsidedict["orf_prot_sequence"] = ""

											query_length = float(((1+sinsidedict['query_end']) - sinsidedict['query_start']) / 3)
											max_ident = float(sinsidedict['max-identities'])
											sinsidedict["perc_identity"] = format(((max_ident*100) / query_length), '.2f')

										elif inType == 'protein':
											sinsidedict["query_start"] = hsp.query_start
											sinsidedict["query_end"] = hsp.query_start + realQueryLength
											sinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
											sinsidedict["orf_prot_sequence"] = orf_protein_sequence
												
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
										slinsidedict["model_type_id"] = modelTypeID
										
										slinsidedict["pass_evalue"] = passevalue
										slinsidedict["model_id"] = modelID
										slinsidedict["ARO_accession"] = json_data[modelID]["ARO_accession"]
										slinsidedict["ARO_name"] = json_data[modelID]["ARO_name"]
										slinsidedict["ARO_category"] = json_data[modelID]["ARO_category"]
										slinsidedict["evalue"] = hsp.expect
										slinsidedict["bit-score"] = hsp.bits
										slinsidedict["max-identities"] = hsp.identities
										if checkKeyExisted("NCBI_taxonomy", json_data[modelID]["model_sequences"]["sequence"][seqinModel]):
											slinsidedict["cvterm_id"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]

										slinsidedict["query"] = hsp.query.encode('ascii','replace')
										slinsidedict["match"] = hsp.match.encode('ascii','replace')
										slinsidedict["sequenceFromDB"] = hsp.sbjct.encode('ascii','replace')
										slinsidedict["SequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"]
										slinsidedict["dnaSequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["dna_sequence"]["sequence"]


										if inType == 'contig':
											if orf == "genemark":
												slinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
												slinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
												slinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
												slinsidedict["orf_start"] = findnthbar(orfInfo, 1)
												slinsidedict["orf_end"] = findnthbar(orfInfo, 2)
												slinsidedict["orf_From"] = orffrom

												if orfInfo[:orfInfo.index('|')] in predicted_genes_dict:
													slinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('|')]]
													slinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
												else:
													slinsidedict["orf_dna_sequence"] = ""
													slinsidedict["orf_prot_sequence"] = ""
											else:
												slinsidedict["query_start"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3
												slinsidedict["query_end"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
												slinsidedict["orf_strand"] = findnthbar2(orfInfo, 3)
												slinsidedict["orf_start"] = findnthbar2(orfInfo, 1)
												slinsidedict["orf_end"] = findnthbar2(orfInfo, 2)
												slinsidedict["orf_From"] = findnthbar2(orfInfo, 0)

												if orfInfo[:orfInfo.index('#')] in predicted_genes_dict:
													slinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('#')]]
													slinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
												else:
													slinsidedict["orf_dna_sequence"] = ""
													slinsidedict["orf_prot_sequence"] = ""

											query_length = float(((1+slinsidedict['query_end']) - slinsidedict['query_start']) / 3)
											max_ident = float(slinsidedict['max-identities'])
											slinsidedict["perc_identity"] = format(((max_ident*100) / query_length), '.2f')

										elif inType == 'protein':
											slinsidedict["query_start"] = hsp.query_start
											slinsidedict["query_end"] = hsp.query_start + realQueryLength
											slinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
											slinsidedict["orf_prot_sequence"] = orf_protein_sequence
												
										loose[hitid + "|hsp_num:" + str(init)] = slinsidedict
										init += 1
		

				elif modelTypeID == 40292:
					init = 0
					passevalue = findnthbar(alignTitle, 1)

					for hsp in alignment.hsps:

						querySeq = hsp.query.replace('-', '')
						realQueryLength = len(querySeq)

						card_sequence = str(json_data[modelID]["model_sequences"]["sequence"][seqinModel]["protein_sequence"]["sequence"])
						orf_protein_sequence = ""
						
						if orf == "genemark":
							if predicted_genes_dict:
								orf_protein_sequence = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
							if submitted_proteins_dict:
								orf_protein_sequence = str(submitted_proteins_dict[orfInfo.split("#")[0]])
						else:
							if predicted_genes_dict:
								if orfInfo in predicted_genes_dict.keys():
									orf_protein_sequence = str(Seq(predicted_genes_dict[orfInfo], generic_dna).translate(table=11)).strip("*")
								else:
									orf_protein_sequence = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
							if submitted_proteins_dict:
								orf_protein_sequence = str(submitted_proteins_dict[orfInfo.split(" ")[0]])														

						if card_sequence == orf_protein_sequence:
							ppinsidedict = {}
							ppinsidedict["type_match"] = "Perfect"
							ppinsidedict["model_id"] = modelID
							ppinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
							ppinsidedict["orf_start"] = findnthbar(orfInfo, 1)
							ppinsidedict["orf_end"] = findnthbar(orfInfo, 2)
							ppinsidedict["orf_From"] = orffrom

							ppinsidedict["model_name"] = json_data[modelID]["model_name"]
							ppinsidedict["model_type"] = json_data[modelID]["model_type"]
							ppinsidedict["model_type_id"] = modelTypeID

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
							ppinsidedict["dnaSequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["dna_sequence"]["sequence"]

							if inType == 'contig':
								if orf == "genemark":
									ppinsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
									ppinsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									ppinsidedict["orf_strand"] = findnthbar(orfInfo, 0)
									ppinsidedict["orf_start"] = findnthbar(orfInfo, 1)
									ppinsidedict["orf_end"] = findnthbar(orfInfo, 2)
									ppinsidedict["orf_From"] = orffrom
									if orfInfo[:orfInfo.index('|')] in predicted_genes_dict:
										ppinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('|')]] 
										ppinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
									else:
										ppinsidedict["orf_dna_sequence"] = ""
										ppinsidedict["orf_prot_sequence"] = ""

								else:
									ppinsidedict["query_start"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3
									ppinsidedict["query_end"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									ppinsidedict["orf_strand"] = findnthbar2(orfInfo, 3)
									ppinsidedict["orf_start"] = findnthbar2(orfInfo, 1)
									ppinsidedict["orf_end"] = findnthbar2(orfInfo, 2)
									ppinsidedict["orf_From"] = findnthbar2(orfInfo, 0)

									if orfInfo[:orfInfo.index('#')] in predicted_genes_dict:
										ppinsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('#')]] 
										ppinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
									else:
										ppinsidedict["orf_dna_sequence"] = ""
										ppinsidedict["orf_prot_sequence"] = ""

								query_length = float(((1+ppinsidedict['query_end']) - ppinsidedict['query_start']) / 3)
								max_ident = float(ppinsidedict['max-identities'])
								ppinsidedict["perc_identity"] = format(((max_ident*100) / query_length), '.2f')
							
							elif inType == 'protein':
								ppinsidedict["query_start"] = hsp.query_start
								ppinsidedict["query_end"] = hsp.query_start + realQueryLength
								ppinsidedict["query_From"] = blast_record.query.encode('ascii','replace')
								ppinsidedict["orf_prot_sequence"] = orf_protein_sequence


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
							insidedict["model_type_id"] = modelTypeID

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
							insidedict["dnaSequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["dna_sequence"]["sequence"]

							if inType == 'contig':
								if orf == "genemark":
									insidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
									insidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									insidedict["orf_strand"] = findnthbar(orfInfo, 0)
									insidedict["orf_start"] = findnthbar(orfInfo, 1)
									insidedict["orf_end"] = findnthbar(orfInfo, 2)
									insidedict["orf_From"] = orffrom
									if orfInfo[:orfInfo.index('|')] in predicted_genes_dict:
										insidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('|')]] 
										insidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
									else:
										insidedict["orf_dna_sequence"] = ""
										insidedict["orf_prot_sequence"] = ""
								else:
									insidedict["query_start"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3
									insidedict["query_end"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									insidedict["orf_strand"] = findnthbar2(orfInfo, 3)
									insidedict["orf_start"] = findnthbar2(orfInfo, 1)
									insidedict["orf_end"] = findnthbar2(orfInfo, 2)
									insidedict["orf_From"] = findnthbar2(orfInfo, 0)
									
									if orfInfo[:orfInfo.index('#')] in predicted_genes_dict:
										insidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('#')]] 
										insidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
									else:
										insidedict["orf_dna_sequence"] = ""
										insidedict["orf_prot_sequence"] = ""									

								query_length = float(((1+insidedict['query_end']) - insidedict['query_start']) / 3)
								max_ident = float(insidedict['max-identities'])
								insidedict["perc_identity"] = format(((max_ident*100) / query_length), '.2f')

							elif inType == 'protein':
								insidedict["query_start"] = hsp.query_start
								insidedict["query_end"] = hsp.query_start + realQueryLength
								insidedict["query_From"] = blast_record.query.encode('ascii','replace')
								insidedict["orf_prot_sequence"] = orf_protein_sequence

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
							linsidedict["model_type_id"] = modelTypeID

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
							linsidedict["dnaSequenceFromBroadStreet"] = json_data[modelID]["model_sequences"]["sequence"][seqinModel]["dna_sequence"]["sequence"]

							if inType == 'contig':
								if orf == "genemark":
									linsidedict["query_start"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3
									linsidedict["query_end"] = findnthbar(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									linsidedict["orf_strand"] = findnthbar(orfInfo, 0)
									linsidedict["orf_start"] = findnthbar(orfInfo, 1)
									linsidedict["orf_end"] = findnthbar(orfInfo, 2)
									linsidedict["orf_From"] = orffrom

									if orfInfo[:orfInfo.index('|')] in predicted_genes_dict:
										linsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('|')]]
										linsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('|')]], generic_dna).translate(table=11)).strip("*")
									else:
										linsidedict["orf_dna_sequence"] = ""
										linsidedict["orf_prot_sequence"] = ""

								else:
									linsidedict["query_start"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3
									linsidedict["query_end"] = findnthbar2(orfInfo, 1) + (hsp.query_start - 1)*3 + realQueryLength*3 - 1
									linsidedict["orf_strand"] = findnthbar2(orfInfo, 3)
									linsidedict["orf_start"] = findnthbar2(orfInfo, 1)
									linsidedict["orf_end"] = findnthbar2(orfInfo, 2)
									linsidedict["orf_From"] = findnthbar2(orfInfo, 0)

									if orfInfo[:orfInfo.index('#')] in predicted_genes_dict:
										linsidedict["orf_dna_sequence"] = predicted_genes_dict[orfInfo[:orfInfo.index('#')]]
										linsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orfInfo[:orfInfo.index('#')]], generic_dna).translate(table=11)).strip("*")
									else:
										linsidedict["orf_dna_sequence"] = ""
										linsidedict["orf_prot_sequence"] = ""

								query_length = float(((1+linsidedict['query_end']) - linsidedict['query_start']) / 3)
								max_ident = float(linsidedict['max-identities'])
								linsidedict["perc_identity"] = format(((max_ident*100) / query_length), '.2f')

							elif inType == 'protein':
								linsidedict["query_start"] = hsp.query_start
								linsidedict["query_end"] = hsp.query_start + realQueryLength
								linsidedict["query_From"] = blast_record.query.encode('ascii','replace')
								linsidedict["orf_prot_sequence"] = orf_protein_sequence
									
							loose[hitid + "|hsp_num:" + str(init)] = linsidedict
							init += 1

			if len(perfect) == 0 and len(strict) == 0:
				if criteria == "no":
					blastResults[blast_record.query.encode('ascii','replace')] = loose
					logging.info("runBlast => hit = " + str(blast_record.query.encode('ascii','replace')) + " => Loose")
				pscore = 1
				
			elif len(perfect) == 0:
				blastResults[blast_record.query.encode('ascii','replace')] = strict
				logging.info("runBlast => hit = " + str(blast_record.query.encode('ascii','replace')) + " => Strict")
				pscore = 2
				
			else:
				blastResults[blast_record.query.encode('ascii','replace')] = perfect
				logging.info("runBlast => hit = " + str(blast_record.query.encode('ascii','replace')) + " => Perfect")
				pscore = 3
		
		blastResults["_metadata"] = {"data_type": data_type}
		pjson = json.dumps(blastResults)

		with open(working_directory+"/"+outputFile+".json", 'w') as f:
			print>>f, pjson
		f.close()

		if inType == 'contig':
			for gene in blastResults:
				logging.debug("runBlast => gene = " + str(gene))
				if not gene == "_metadata":
					for hsp in blastResults[gene]:
						if blastResults[gene][hsp]["orf_end"] < blastResults[gene][hsp]["query_end"]:
							pass							
							#print>>sys.stderr, hsp
							#print>>sys.stderr, blastResults[gene][hsp]["orf_end"]
							#print>>sys.stderr, blastResults[gene][hsp]["query_end"]
							#logging.error("runBlast => hsp => " + str(hsp))
							#logging.error("runBlast => [gene][hsp] => " + str(blastResults[gene][hsp]))
	logging.info("runBlast => done")
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

def decompress(inputfile,ext,out_dir):

    filename = os.path.basename(inputfile)
    name = os.path.splitext(filename)
    newFile = out_dir + "/" + str(name[0])
    clean_files.append(newFile)

    if ext == 'gz':
	    readFq=gzip.open(inputfile)
	    data = readFq.readlines()

	    with open(newFile,'w') as fin:
	    	for eachline in data:
	    		print>>fin, eachline.rstrip() 
	    fin.close()

    return newFile

def validateHeaders(inputfile,orf):
	message = []
	from Bio import SeqIO
	for record in SeqIO.parse(inputfile, 'fasta'):
		if "#" in record.id and orf == "prodigal":
			# Bad header with character #
			message.append("Bad header with character # : >"+ str(record.id))
		elif "|" in record.id and orf == "genemark":
			# Bad header with character | 
			message.append("Bad header with character | : >"+ str(record.id))

	if message:
		#print message
		print>>sys.stderr, message
		exit()


def main(args):
	if args.inputSeq == None:
		print "[error] missing input(s)"
		print "[info] Try: rgi -h"
		logging.error("main => missing input(s)")
		logging.error("main => Try: rgi -h")
		exit()
	if args.verbose.lower() == 'on':
		log_file = str(working_directory) + "/"+os.path.basename(args.inputSeq)+".log"
		print "[info] logs saved to : " + log_file
		logging.basicConfig(filename=log_file, level=logging.DEBUG, format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')
		logging.info('main => start rgi')
	
	inType = args.inType.lower()
	inputSeq = args.inputSeq 
	threads = args.threads
	outputFile = args.output
	criteria = args.criteria.lower()
	clean = args.clean.lower()
	data_type = args.data.lower()
	orf = args.orf.lower()
	alignment_tool = args.aligner.lower()

	bpjson = ""

	# Write each request to a directory based on the input filename with (-results) appeneded
	#output_dir = working_directory + "/" + os.path.basename(inputSeq) + "-results"
	output_dir = working_directory

	#if not os.path.exists(output_dir):
	#	os.makedirs(output_dir)

	# Check if file is compressed
	if inputSeq.endswith('.gz'):
		inputSeq = decompress(inputSeq,'gz',output_dir)

	checkBeforeBlast(inType, inputSeq)
	writeFASTAfromJson()

	if alignment_tool == "blast":
		makeBlastDB(inType, inputSeq)
	else:
		makeDiamondDB()

	# check the heasders
	validateHeaders(inputSeq,orf)
	# run blast
	bpjson = runBlast(inType, inputSeq, threads, outputFile, criteria, data_type, clean, orf, alignment_tool)

	# Generate tab-delimited files
	logging.info('main => Generate tab-delimited ' + outputFile + ".txt")
	convertJsonToTSV.printCSV(outputFile + ".json",outputFile,orf,args.verbose.lower())

	if clean == "yes":
		logging.info('main => clean temporary files')
		removeTemp()
	logging.info('main => rgi success')
	logging.info('----------------------------------------')
	return bpjson

def version():
	software_version = "3.1.2"
	return software_version

def data_version():
	data_version = ""
	with open(data_path+"card.json") as json_file:
		json_data = json.load(json_file)
		for item in json_data.keys():
			if item == "_version":
				data_version = json_data[item]
	return data_version

def run():
	parser = argparse.ArgumentParser(description='Resistance Gene Identifier - Version ' + version())
	parser.add_argument('-t','--input_type', dest="inType", default="contig", help='must be one of contig, orf, protein, read (default: contig)')
	parser.add_argument('-i','--input_sequence',dest="inputSeq",help='input file must be in either FASTA (contig and protein), FASTQ(read) or gzip format! e.g myFile.fasta, myFasta.fasta.gz')
	parser.add_argument('-n', '--num_threads',  dest="threads", default="32", help="Number of threads (CPUs) to use in the BLAST search (default=32)")
	parser.add_argument('-o', '--out_file',  dest="output", default="Report", help="Output JSON file (default=Report)")	
	parser.add_argument('-e', '--exclude_loose',  dest="criteria", default="YES", help="This option is used to include or exclude the loose hits. Options are NO or YES (default=YES for excluding loose hits)")
	parser.add_argument('-c', '--clean',  dest="clean", default="YES", help="This removes temporary files in the results directory after run. Options are NO or YES (default=YES for remove)")
	parser.add_argument('-d', '--data', dest="data", default="NA", help = "Specify a data-type, i.e. wgs, chromosome, plasmid, etc. (default = NA)")
	parser.add_argument('-l', '--verbose', dest="verbose", default="OFF", help = "log progress to file. Options are OFF or ON  (default = OFF for no logging)")
	parser.add_argument('-x', '--orf', dest="orf", default="PRODIGAL", help = "choose between prodigal and MetaGeneMark orf finder. Options are PRODIGAL or GENEMARK  (default = PRODIGAL)")
	parser.add_argument('-a', '--alignment_tool', dest="aligner", default="BLAST", help = "choose between BLAST and DIAMOND. Options are BLAST or DIAMOND  (default = BLAST)")
	parser.add_argument('-sv', '--software_version', action='version', version='rgi-software-'+version(), help = "Prints software number")
	parser.add_argument('-dv', '--data_version', action='version', version='rgi-data-version-'+data_version() ,help = "Prints data version number")
	args = parser.parse_args()
	main(args)	

if __name__ == '__main__':
	run()
	
