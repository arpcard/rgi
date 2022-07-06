from app.Base import BaseModel
from app.settings import *

class Variant(BaseModel):
	"""Class for protein variant searches."""
	def __init__(self, input_type, loose, input_sequence, xml_file, working_directory, local_database=False, include_nudge=False):
		self.input_type = input_type
		self.loose = loose
		self.input_sequence = input_sequence
		self.xml_file = xml_file
		self.output = {}
		self.working_directory = working_directory

		self.local_database = local_database
		self.data = data_path

		self.include_nudge = include_nudge

		if self.local_database:
			self.db = LOCAL_DATABASE
			self.data = LOCAL_DATABASE

	def __repr__(self):
		"""Returns Variant class full object."""
		return "Variant({}".format(self.__dict__)

	def run(self):
		blastResults = {}
		predicted_genes_dict = {}
		predicted_genes_dict_protein = {}
		submitted_proteins_dict = {}
		orf=0

		if self.input_type == "contig":
			predicted_genes_dict = self.get_orf_dna_sequence(self.input_sequence,self.input_type)
			predicted_genes_dict_protein = self.get_orf_protein_sequence(self.input_sequence,self.input_type)

		if self.input_type == "protein":
			submitted_proteins_dict = (self.get_submitted_protein_sequence(self.input_sequence))

		with open(os.path.join(self.data,"card.json")) as json_file:
			json_data = json.load(json_file)

		with open(self.xml_file, 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)
			for blast_record in blast_records:
				perfect = {}
				strict = {}
				loose = {}
				for alignment in blast_record.alignments:
					align_title = alignment.title
					orf_info = blast_record.query.encode('ascii','replace')

					c = 0
					barc = 0
					for eachc in orf_info:
						if barc >= 6:
							break
						elif eachc == '|':
							barc += 1
							c += 1
						else:
							c += 1
					orf_from = orf_info[c:]
					
					model_type_id = self.extract_nth_bar(align_title, 0)
					# logger.info("model_type_id: {} ".format(model_type_id))
					space_pos = align_title.index(' ')
					hit_id = align_title[0:space_pos]
					hit_id = hit_id.encode('ascii','replace')
					model_descrpt =align_title[align_title.index(' ')+1:]
					underscore_in_MD = model_descrpt.index('_')
					model_id = model_descrpt[0:underscore_in_MD]
					seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
					pass_value = self.extract_nth_bar(alignment.title, 1)
					# logger.info("pass_value: {}".format(pass_value))

					if model_type_id == 40293:
						try:
							true_pass_evalue = float(pass_value)
						except ValueError:
							true_pass_evalue = float(pass_value[0:pass_value.find(' ')])

						# logger.info("mutation | model_type_id = " + str(align_title))				
						init = 0
						evalue_snp = self.extract_nth_bar(align_title, 2)
						snpl = []
						snp_dict_list = []
						temp = ""
						evalue_snp_dec = evalue_snp
						snpl = evalue_snp_dec.split(',')
						
						for each_snp in snpl:
							snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})

						for hsp in alignment.hsps:
							query_seq =  hsp.query.replace('-', '')
							real_query_length = len(query_seq) 
							sbjct_seq = hsp.sbjct.replace('-', '') 
							real_sbjct_length = len(sbjct_seq) 

							for eachs in snp_dict_list:
								pos = eachs["position"]
								ori = eachs["original"]
								chan = eachs["change"]

								if hsp.sbjct_start < pos and (hsp.sbjct_start + real_sbjct_length) > pos:
									orf_protein_sequence = ""

									# if predicted_genes_dict:
									# 	if orf_info.strip() in predicted_genes_dict.keys():
									# 		orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()], generic_dna).translate(table=11)).strip("*")
									# 	else:
									# 		orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()[:orf_info.decode().index(' # ')]], generic_dna).translate(table=11)).strip("*")

									if predicted_genes_dict_protein:
										if orf_info.strip() in predicted_genes_dict_protein.keys():
											orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()].strip("*")
										else:
											orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")

									if submitted_proteins_dict:
										orf_protein_sequence = str(submitted_proteins_dict[orf_info.decode().split(" ")[0]])
										
									# logger.info("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
									# 			self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
									# 			+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

									# Report ONLY if the SNPs are present
									qry = int(pos) - hsp.sbjct_start + self.find_num_dash(hsp.sbjct, (int(pos) - hsp.sbjct_start))
									sbj = int(pos) - hsp.sbjct_start + self.find_num_dash(hsp.sbjct, (int(pos) - hsp.sbjct_start))

									if hsp.query[qry] == chan:
										query_snps = {}	
										# logger.debug("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
										# 		self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
										# 		+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

										# get position of mutation in the query sequence
										d = int(pos) - hsp.sbjct_start - self.find_num_dash(hsp.query, (int(pos) - hsp.sbjct_start))
										# print("hsp.sbjct_start: ", hsp.sbjct_start)
										query_snps = {"original": ori, "change": chan ,"position": d+1}
										# logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

										try:
											if float(hsp.bits) >= float(true_pass_evalue):		
												sinsidedict = {}
												sinsidedict["type_match"] = "Strict"
												sinsidedict["snp"] = eachs
												sinsidedict["query_snp"] = query_snps
												sinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
												sinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)
												sinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
												sinsidedict["orf_from"] = orf_from.decode()
												sinsidedict["model_name"] = json_data[model_id]["model_name"]
												sinsidedict["model_type"] = json_data[model_id]["model_type"]
												sinsidedict["model_type_id"] = model_type_id
												sinsidedict["model_id"] = model_id
												sinsidedict["pass_evalue"] = "n/a"
												sinsidedict["pass_bitscore"] = pass_value
												sinsidedict["ARO_accession"] = json_data[model_id]["ARO_accession"]
												sinsidedict["ARO_name"] = json_data[model_id]["ARO_name"]
												sinsidedict["ARO_category"] = json_data[model_id]["ARO_category"]
												sinsidedict["evalue"] = hsp.expect
												sinsidedict["max_identities"] = hsp.identities
												sinsidedict["bit_score"] = hsp.bits
												sinsidedict["cvterm_id"]  = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
												sinsidedict["query"] = hsp.query
												sinsidedict["match"] = hsp.match
												sinsidedict["sequence_from_db"] = hsp.sbjct
												sinsidedict["sequence_from_broadstreet"]	= json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												sinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													sinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													sinsidedict["partial"] = "0"
												
												if self.input_type == 'contig':
													sinsidedict["query_start"] = self.extract_nth_hash(orf_info.decode(), 1) + (hsp.query_start - 1)*3
													sinsidedict["query_end"] = self.extract_nth_hash(orf_info.decode(), 1) + (hsp.query_start - 1)*3 + real_query_length*3 - 1
													sinsidedict["orf_strand"] = self.extract_nth_hash(orf_info.decode(), 3)
													sinsidedict["orf_start"] = self.extract_nth_hash(orf_info.decode(), 1)
													sinsidedict["orf_end"] = self.extract_nth_hash(orf_info.decode(), 2)
													sinsidedict["orf_from"] = self.extract_nth_hash(orf_info.decode(), 0)

													if orf_info.decode().split(' # ')[0] in predicted_genes_dict:
														sinsidedict["orf_dna_sequence"] = predicted_genes_dict[orf_info.decode().split(' # ')[0]] 
														# sinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orf_info.decode().split(' # ')[0]], generic_dna).translate(table=11)).strip("*")
														sinsidedict["orf_prot_sequence"] = orf_protein_sequence
													else:
														sinsidedict["orf_dna_sequence"] = ""
														sinsidedict["orf_prot_sequence"] = ""	


												elif self.input_type == 'protein':
													sinsidedict["query_start"] = hsp.query_start
													sinsidedict["query_end"] = hsp.query_start + real_query_length
													sinsidedict["query_from"] = blast_record.query
													sinsidedict["orf_prot_sequence"] = orf_protein_sequence

												elif self.input_type == 'read':
													pass

												sinsidedict["perc_identity"] = float(format(float(sinsidedict["max_identities"]*100) / len(sinsidedict["query"]),'.2f'))

												strict["{}|hsp_num:{}".format(hit_id.decode(),init)] = sinsidedict
												init += 1

											else:
												slinsidedict = {}
												slinsidedict["type_match"] = "Loose"
												slinsidedict["snp"] = eachs
												slinsidedict["query_snp"] = query_snps
												slinsidedict["orf_strand"] = self.extract_nth_bar(orf_info.decode(), 0)
												slinsidedict["orf_start"] = self.extract_nth_bar(orf_info.decode(), 1)				
												slinsidedict["orf_end"] = self.extract_nth_bar(orf_info.decode(), 2)
												slinsidedict["orf_from"] = orf_from.decode()
												slinsidedict["model_name"] = json_data[model_id]["model_name"]
												slinsidedict["model_type"] = json_data[model_id]["model_type"]
												slinsidedict["model_type_id"] = model_type_id
												slinsidedict["pass_evalue"] = "n/a"
												slinsidedict["pass_bitscore"] = pass_value
												slinsidedict["model_id"] = model_id
												slinsidedict["ARO_accession"] = json_data[model_id]["ARO_accession"]
												slinsidedict["ARO_name"] = json_data[model_id]["ARO_name"]
												slinsidedict["ARO_category"] = json_data[model_id]["ARO_category"]
												slinsidedict["evalue"] = hsp.expect
												slinsidedict["bit_score"] = hsp.bits
												slinsidedict["max_identities"] = hsp.identities
												slinsidedict["cvterm_id"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["NCBI_taxonomy"]["NCBI_taxonomy_cvterm_id"]
												slinsidedict["query"] = hsp.query
												slinsidedict["match"] = hsp.match
												slinsidedict["sequence_from_db"] = hsp.sbjct
												slinsidedict["sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												slinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													slinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													slinsidedict["partial"] = "0"

												if self.input_type == 'contig':
													slinsidedict["query_start"] = self.extract_nth_hash(orf_info.decode(), 1) + (hsp.query_start - 1)*3
													slinsidedict["query_end"] = self.extract_nth_hash(orf_info.decode(), 1) + (hsp.query_start - 1)*3 + real_query_length*3 - 1
													slinsidedict["orf_strand"] = self.extract_nth_hash(orf_info.decode(), 3)
													slinsidedict["orf_start"] = self.extract_nth_hash(orf_info.decode(), 1)
													slinsidedict["orf_end"] = self.extract_nth_hash(orf_info.decode(), 2)
													slinsidedict["orf_from"] = self.extract_nth_hash(orf_info.decode(), 0)

													if orf_info.decode().split(' # ')[0] in predicted_genes_dict:
														slinsidedict["orf_dna_sequence"] = predicted_genes_dict[orf_info.decode().split(' # ')[0]] 
														# slinsidedict["orf_prot_sequence"] = str(Seq(predicted_genes_dict[orf_info.decode().split(' # ')[0]], generic_dna).translate(table=11)).strip("*")
														slinsidedict["orf_prot_sequence"] = orf_protein_sequence
													else:
														slinsidedict["orf_dna_sequence"] = ""
														slinsidedict["orf_prot_sequence"] = ""	

												elif self.input_type == 'protein':
													slinsidedict["query_start"] = hsp.query_start
													slinsidedict["query_end"] = hsp.query_start + real_query_length
													slinsidedict["query_from"] = blast_record.query
													slinsidedict["orf_prot_sequence"] = orf_protein_sequence

												elif self.input_type == 'read':
													pass

												slinsidedict["perc_identity"] = float(format(float(slinsidedict["max_identities"]*100) / len(slinsidedict["query"]),'.2f'))
												loose["{}|hsp_num:{}".format(hit_id.decode(),init)] = slinsidedict

												init += 1
										except Exception as e:
											logger.warning("Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
											logger.warning("{} ---> hsp.bits: {} {} ? {}".format(json_data[model_id]["model_name"],hsp.bits,type(hsp.bits), type(true_pass_evalue)))

				blastResults = self.results(blastResults, blast_record.query, perfect, strict , loose, self.include_nudge)
						
			return blastResults
	