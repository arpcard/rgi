import sys
import os
import filepaths
from Bio import SeqIO
from Bio.Seq import Seq
import json

'''
Bacteria (70S)  
        LSU 50S
                5S      RF00001
                23S     SILVA-LSU-Bac
        SSU 30S
                16S     RF00177
'''

'''
Step 1: get rRNA coordinates using barrnap
Step 2: extract rRNA sequences
Step 3: blast and report for mutations in CARD
'''

script_path = filepaths.determine_path()
working_directory = os.getcwd()

path = script_path+"/_db/"
data_path = script_path+"/_data/"
tmp = script_path+"/_tmp/"

def rRNA(fq_path,start,stop,strand):
	try:
		sequence = SeqIO.read(fq_path, "fasta")
		return str(sequence.seq[int(start):int(stop)])
	except Exception as e:
		print "Error: ",e, ".Please input only one record / sequence"
		return ""

def rRNA_list(fq_path,file_name):
	list_all = {}
	list_all["5s"] = []
	list_all["16s"] = []
	list_all["23s"] = []
	
	# get rRNA
	os.system('barrnap --reject 1e-9 --quiet --kingdom bac '+fq_path+' > '+file_name+'.gff3') # --incseq

	# scan file and get sequence and coordinates
	with open(working_directory+'/'+file_name+'.gff3', 'r') as f:
		data = f.readlines()
		for eachline in data:
			if eachline == '':
				pass
			elif eachline[0:13] == '##gff-version':
				pass
			elif 'Name=5S_rRNA;product=5S ribosomal RNA' in eachline:
				list_5s = eachline.split("\t")
				list_all["5s"].append({"start": list_5s[3], "stop": list_5s[4], "strand": str(list_5s[6]), "sequence": rRNA(fq_path,list_5s[3],list_5s[4],list_5s[6])})
			elif 'Name=16S_rRNA;product=16S ribosomal RNA' in eachline:
				list_16s = eachline.split("\t")
				list_all["16s"].append({"start": list_16s[3], "stop": list_16s[4], "strand": str(list_16s[6]), "sequence": rRNA(fq_path,list_16s[3],list_16s[4],list_16s[6])})
			elif 'Name=23S_rRNA;product=23S ribosomal RNA' in eachline:
				list_23s = eachline.split("\t")
				list_all["23s"].append({"start": list_23s[3], "stop": list_23s[4], "strand": str(list_23s[6]), "sequence": rRNA(fq_path,list_23s[3],list_23s[4],list_23s[6])})

	pjson = json.dumps(list_all)

	with open(working_directory+"/"+file_name+".rRNA.json", 'w') as f:
		print>>f, pjson
	f.close()

def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone

def writeFASTAfromJson():
	noSeqList = []
	if os.path.isfile(path+"nucleotidedb.fsa") == False:
		with open(data_path+"card.json") as json_file:
			json_data = json.load(json_file)
			with open (path+'nucleotidedb.fsa', 'w') as wp:
				for item in json_data:
					if item.isdigit():
						if json_data[item]["model_type_id"] == "40295":
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
									print>>wp, ('>' + item + '_' + seqkey + " | model_type_id: 40295" + " | pass_evalue: " + str(pass_eval) + " | SNP: " + snpList)
									print>>wp, (json_data[item]["model_sequences"]["sequence"][seqkey]["dna_sequence"]["sequence"])
							else:
								 noSeqList.append(item)
					
		json_file.close()
	print noSeqList

def runDiamond(inType, inputSeq, threads, outputFile):

	cfilter = "	--index-chunks 1 --block-size 1 --quiet"

	if inType == 'contig':
		#logging.info("runDiamond => blastx")
		#logging.info('runDiamond => diamond blastx --in '+path+'nucleotidedb.fsa --db '+path+'nucleotidedb.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles --more-sensitive --evalue 10'+cfilter)
		os.system('diamond blastx --in '+path+'nucleotidedb.fsa --db '+path+'nucleotidedb.db'+' --query '+inputSeq+' --daa '+inputSeq+' --threads '+threads+' --salltitles --more-sensitive --evalue 10'+cfilter)
	else:
		print "Error"
    
	#logging.info('runDiamond => diamond view --outfmt xml --daa '+inputSeq+'.daa --out '+outputFile + cfilter)
	os.system('diamond view --outfmt xml --daa '+inputSeq+'.daa --out '+outputFile + cfilter )
	#logging.info('runDiamond => diamond view --daa '+inputSeq+'.daa --out '+outputFile+".txt" + cfilter)
	os.system('diamond view --daa '+inputSeq+'.daa --out '+outputFile+".txt" + cfilter )

	#clean_files.append(inputSeq+'.daa')
	#clean_files.append(outputFile+'.txt')
def runBlast(inputSeq,threads):
	file_name = os.path.basename(inputSeq)
	os.system("blastn -perc_identity 0.10 -out blastnRes.xml -outfmt 5 -query " + inputSeq + " -db "+path+"dna.db")


def makeDiamondDB():
	if os.path.isfile(path+"nucleotidedb.fsa") == True and os.path.exists(path+"nucleotidedb.fsa") == True  and os.path.exists(path+"nucleotidedb.db.dmnd") == True:
		print "diamond DB exists"
	else:
		print "create diamond DB."
		os.system('diamond makedb --in '+path+'nucleotidedb.fsa --db '+path+'nucleotidedb.db')

def makeBlastDB():
	if os.path.isfile(path+"nucleotidedb.fsa") == True and os.path.exists(path+"nucleotidedb.fsa") == True  and os.path.exists(path+"dna.db.phr") == True and os.path.exists(path+"dna.db.pin") == True and os.path.exists(path+"dna.db.psq") == True :
		print "blast DB exists"
	else:
		print "create blast DB."
		os.system('makeblastdb -in '+path+'nucleotidedb.fsa -dbtype nucl -out '+path+'dna.db')
	#os.system('makeblastdb -in proteindb.fsa -dbtype prot -out protein.db')

def main(argvfile):
	exit("In development...")
	#file_name = os.path.basename(argvfile)
	# get rRNA sequences
	#rRNA_list(argvfile,file_name)
	#rRNA(argvfile,0,1,'+')
	#writeFASTAfromJson()
	#makeDiamondDB()
	#runDiamond('contig', argvfile, '8', file_name)
	#makeBlastDB()
	#runBlast(argvfile,8)

if __name__ == "__main__":
    main(sys.argv[1])