from app.settings import *

class Output(object):
	"""Class to format RGI output."""
	def __init__(self,input_file):
		"""Creates Output object for formatting outputs."""
		self.input_file = input_file
		f_name = os.path.basename(input_file)
		self.output_file = os.path.join(working_directory,"{}.tab.txt".format(f_name))
		logger.debug('Created Output object')
		logger.info(repr(self))
	def __repr__(self):
		"""Returns Output class full object."""
		return "Output({}".format(self.__dict__)

	def run(self):
		logger.info("Create tab-delimited file.")
		self.print_csv()

	def checkKeyExisted(self,key, my_dict):
		try:
			nonNone = my_dict[key] is not None
		except KeyError:
			nonNone = False
		return nonNone

	#output the information particular field from alignment.Title by splicing it by '#'
	def findnthbar2(self,bunchstr, n):
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

	def print_csv(self):
		if os.path.isfile(self.input_file) == False:
			logger.error("missing input JSON file. ({})".format(self.input_file))
			exit()

		try:
			with open(self.input_file, 'r') as f:
				data = json.load(f)
		
		except ValueError:
			logger.error("invalid JSON string.")
			exit()

		with open(self.output_file, "w") as af:
			writer = csv.writer(af, delimiter='\t', dialect='excel')

			writer.writerow(["ORF_ID", "CONTIG", "START", "STOP", "ORIENTATION", \
							 "CUT_OFF", "PASS_EVALUE", "Best_Hit_evalue", "Best_Hit_ARO", \
							 "Best_Identities", "ARO", "ARO_name", "Model_type", "SNP", \
							 "Best_Hit_ARO_category", "ARO_category", "PASS_bitscore", \
							 "Best_Hit_bitscore", "bit_score","Predicted_DNA","Predicted_Protein",\
							 "CARD_Protein_Sequence","LABEL","ID","Model_id"])
			
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
						if self.checkKeyExisted("ARO_category", data[item][it]):
							for aroctkey in data[item][it]["ARO_category"]:
								#cgList.append(str(data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace')))
								cgList.append(str(data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace'),'utf-8'))

						if data[item][it]["model_type_id"] == 40293:
							temp = data[item][it]["SNP"]["original"] + str(data[item][it]["SNP"]["position"]) + data[item][it]["SNP"]["change"]
							snpList.append(temp)
						elif data[item][it]["model_type_id"] == 40292:
							snpList.append("n/a")

						AROlist.append(str(data[item][it]["ARO_accession"].encode('ascii','replace'), 'utf-8'))
						AROnameList.append(str(data[item][it]["ARO_name"].encode('ascii','replace'), 'utf-8'))
						bitScoreList.append(float(data[item][it]["bit_score"]))
						pass_evalue = "n/a"
						pass_bitscore = float(data[item][it]["pass_bitscore"])
						AROcatList.append(cgList)
						typeList.append(data[item][it]["model_type"])
						cutoffList.append(data[item][it]["type_match"])
						identityList.append(float(data[item][it]["perc_identity"]))
						bestAROcategory = []

						# sort results by max bit-score and maximum percent identity
						if startCompare:
							if maxscore < float(data[item][it]["bit-score"]) and maxpercent < float(data[item][it]["perc_identity"]):
								minevalue = float(data[item][it]["evalue"])
								maxscore = float(data[item][it]["bit_score"])
								maxpercent = float(data[item][it]["perc_identity"])
								minARO = data[item][it]["ARO_name"]
								type_match = data[item][it]["type_match"]
								topModel = data[item][it]["model_id"]
								SequenceFromBroadStreet = data[item][it]["sequence_from_broadstreet"]

								if "orf_prot_sequence" in data[item][it]:
									predictedProtein = data[item][it]["orf_prot_sequence"]
								if "orf_dna_sequence" in data[item][it]:
									predictedDNA = data[item][it]["orf_dna_sequence"]

								if self.checkKeyExisted("ARO_category", data[item][it]):
									for key in data[item][it]["ARO_category"]:
										bestAROcategory.append(str(data[item][it]["ARO_category"][key]["category_aro_name"].encode('ascii','replace'), 'utf-8'))
									bestAROcategorydict[str(minARO)+"|"+str(minevalue)] = bestAROcategory

								if "hsp_num:" in it:
									hitID = it
							
						else:
							startCompare = True
							minevalue = data[item][it]["evalue"]
							maxscore = int(data[item][it]["bit_score"])
							maxpercent = float(data[item][it]["perc_identity"])
							minARO = data[item][it]["ARO_name"]
							type_match = data[item][it]["type_match"]
							topModel = data[item][it]["model_id"]
							SequenceFromBroadStreet = data[item][it]["sequence_from_broadstreet"]

							if "orf_prot_sequence" in data[item][it]:
								predictedProtein = data[item][it]["orf_prot_sequence"]
							if "orf_dna_sequence" in data[item][it]:
									predictedDNA = data[item][it]["orf_dna_sequence"]

							if self.checkKeyExisted("ARO_category", data[item][it]):
								for key in data[item][it]["ARO_category"]:
									bestAROcategory.append(str(data[item][it]["ARO_category"][key]["category_aro_name"].encode('ascii','replace'), 'utf-8'))
								bestAROcategorydict[str(minARO)+"|"+str(minevalue)] = bestAROcategory

							if "hsp_num:" in it:
								hitID = it
							print("AMOS", minARO, "-> ", type_match)

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
						if self.findnthbar2(item, 1) == "":
							writer.writerow([item, 
								"", 
								"", 
								"", 
								"", 
								#', '.join(list(clist)),
								type_match,
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
							writer.writerow([self.findnthbar2(item, 0),
								self.findnthbar2(item, 4).strip(" "), 
					        	int(self.findnthbar2(item, 1))-1, 
					        	int(self.findnthbar2(item, 2))-1, 
					        	self.findnthbar2(item, 3), 
					        	#', '.join(list(clist)), 
					        	type_match,
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


# obj = Output("./_tests/protein.fasta.blast.output.json")
# print(repr(obj))
# obj.run()

# obj = Output("./_tests/nucleotide.fa.diamond.output.json")
# print(repr(obj))
# obj.run()
