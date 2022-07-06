from app.Base import BaseModel
from app.settings import *

class Rrna(BaseModel):
	"""Class for ribosomal RNA searches."""
	def __init__(self, input_file, output_file, db, xml, loose, local_database=False, include_nudge=False, num_threads=32):
		self.input_file = input_file
		self.output_file = output_file
		self.db = db
		self.xml_file = xml
		self.loose = loose

		self.local_database = local_database
		self.data = data_path

		self.include_nudge = include_nudge
		self.num_threads = num_threads

		if self.local_database:
			# self.db = LOCAL_DATABASE
			self.data = LOCAL_DATABASE

	def __repr__(self):
		"""Returns Ribosomal RNA class full object."""
		return "Rrna({}".format(self.__dict__)

	def sequence_orientation(self, end, start):
		if end > start:
			return "+"
		else: 
			return "-"

	def run(self):
		blastResults = {}
		with open(os.path.join(self.data,"card.json")) as json_file:
			json_data = json.load(json_file)

		with open(os.path.join(self.xml_file), 'r') as result_handle:
			blast_records = NCBIXML.parse(result_handle)

			for blast_record in blast_records:
				perfect = {}
				strict = {}
				loose = {}

				for alignment in blast_record.alignments:
					align_title = alignment.title
					orf_info = blast_record.query

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
					orf_info_str = str(orf_info).split(" | ")
					model_type_id = int(orf_info_str[1].split(":")[1].strip())
					# logger.debug("model_type_id: {} ".format(model_type_id))
				
					space_pos = align_title.index(' ')

					hit_id = align_title[0:space_pos]
					hit_id = hit_id.encode('ascii','replace')

					model_descrpt = "?model_descrpt?"
					model_info = orf_info_str[0].split("_")
					model_id = model_info[0].strip()
					seq_in_model = model_info[1].strip()
					pass_value = orf_info_str[2].split(":")[1].strip()

					# logger.debug("model_id: {}".format(model_id))
					# logger.debug("pass_value: {}".format(pass_value))
					
					if model_type_id == 40295:
						true_pass_evalue = float(pass_value)

						init = 0
						evalue_snp = orf_info_str[3].split(":")[1].strip()
						snpl = []
						snp_dict_list = []
						temp = ""
						snpl = evalue_snp.split(',')

						for each_snp in snpl:
							snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})

						for hsp in alignment.hsps:
							query_seq =  hsp.query.replace('-', '')
							real_query_length = len(query_seq) 
							sbjct_seq = hsp.sbjct.replace('-', '') 
							real_sbjct_length = len(sbjct_seq) 
							strand = self.sequence_orientation(hsp.sbjct_end, hsp.sbjct_start)

							for eachs in snp_dict_list:
								pos = eachs["position"]
								ori = eachs["original"]
								chan = eachs["change"]

								if hsp.query_start < pos and (hsp.query_start + real_query_length) > pos:
									# Report ONLY if the SNPs are present								
									qry = int(pos) - hsp.query_start + self.find_num_dash(hsp.query, (int(pos) - hsp.query_start))
									sbj = int(pos) - hsp.query_start + self.find_num_dash(hsp.query, (int(pos) - hsp.query_start))

									if hsp.sbjct[sbj] == chan:
										query_snps = {}
										logger.info("hsp.query_start: {}".format(hsp.query_start))
										logger.info("hsp.query_end: {}".format(hsp.query_end))
										logger.info("hsp.sbjct_start: {}".format(hsp.sbjct_start))
										logger.info("hsp.sbjct_end: {}".format(hsp.sbjct_end))
										d = 0
										if strand == "+":
											d = int(pos) - hsp.query_start - self.find_num_dash(hsp.sbjct, (int(pos) - hsp.query_start))
											# d = -1*(hsp.query_start - pos)
											query_snps = {"original": hsp.query[d], "change": hsp.sbjct[d], "position": (d + 1)}
										else:
											d = int(pos) - hsp.query_start - self.find_num_dash(hsp.sbjct, (int(pos) - hsp.query_start))
											# d = -1*(hsp.query_start - pos)
											query_snps = {"original": hsp.query[d], "change": hsp.sbjct[d], "position": (d + 1)}

										logger.info("position in the query (3'->5') : {}".format(d))

										# logger.info("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))
										try:
											if float(hsp.bits) >= float(true_pass_evalue):
												sinsidedict = {}
												sinsidedict["type_match"] = "Strict"
												sinsidedict["snp"] = eachs
												sinsidedict["query_snp"] = query_snps
												sinsidedict["orf_strand"] = strand
												sinsidedict["orf_start"] = hsp.sbjct_start
												sinsidedict["orf_end"] = hsp.sbjct_end

												sinsidedict["_orf_strand"] = self.extract_nth_bar(orf_info, 0)
												sinsidedict["_orf_start"] = self.extract_nth_bar(orf_info, 1)
												sinsidedict["_orf_end"] = self.extract_nth_bar(orf_info, 2)

												sinsidedict["orf_from"] = alignment.hit_def
												sinsidedict["strand"] = strand
												sinsidedict["hit_def"] = alignment.hit_def
												sinsidedict["sbjct_start"] = hsp.sbjct_start
												sinsidedict["sbjct_end"] = hsp.sbjct_end
												sinsidedict["query_start"] = hsp.query_start
												sinsidedict["query_end"] = hsp.query_end
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
												sinsidedict["sbjct"] = hsp.sbjct
												sinsidedict["sequence_from_db"] = ""
												sinsidedict["orf_dna_sequence"] = hsp.sbjct
												sinsidedict["orf_prot_sequence"] = ""

												sinsidedict["sequence_from_broadstreet"]	= json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												sinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													sinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													sinsidedict["partial"] = "0"
												sinsidedict["perc_identity"] = float(format(float(sinsidedict["max_identities"]*100) / len(sinsidedict["query"]),'.2f'))

												strict["{}|hsp_num:{}".format(hit_id.decode(),init)] = sinsidedict
												init += 1

											else:
												slinsidedict = {}
												slinsidedict["type_match"] = "Loose"
												slinsidedict["snp"] = eachs
												slinsidedict["query_snp"] = query_snps
												slinsidedict["orf_strand"] = strand
												slinsidedict["orf_start"] = hsp.sbjct_start
												slinsidedict["orf_end"] = hsp.sbjct_end

												slinsidedict["_orf_strand"] = self.extract_nth_bar(orf_info, 0)
												slinsidedict["_orf_start"] = self.extract_nth_bar(orf_info, 1)
												slinsidedict["_orf_end"] = self.extract_nth_bar(orf_info, 2)

												slinsidedict["orf_from"] = alignment.hit_def
												slinsidedict["strand"] = strand
												slinsidedict["hit_def"] = alignment.hit_def
												slinsidedict["sbjct_start"] = hsp.sbjct_start
												slinsidedict["sbjct_end"] = hsp.sbjct_end
												slinsidedict["query_start"] = hsp.query_start
												slinsidedict["query_end"] = hsp.query_end
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
												slinsidedict["sequence_from_db"] = ""
												slinsidedict["sbjct"] = hsp.sbjct
												slinsidedict["orf_dna_sequence"] = hsp.sbjct
												slinsidedict["orf_prot_sequence"] = ""

												slinsidedict["sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["protein_sequence"]["sequence"]
												slinsidedict["dna_sequence_from_broadstreet"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["sequence"]
												if "partial" in json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"].keys():
													slinsidedict["partial"] = json_data[model_id]["model_sequences"]["sequence"][seq_in_model]["dna_sequence"]["partial"]
												else:
													slinsidedict["partial"] = "0"
												slinsidedict["perc_identity"] = float(format(float(slinsidedict["max_identities"]*100) / len(slinsidedict["query"]),'.2f'))

												loose["{}|hsp_num:{}".format(hit_id.decode(),init)] = slinsidedict
												init += 1
										except Exception as e:
											logger.warning("Exception : {} -> {} -> Model({})".format(type(e), e, model_id))
											logger.warning("{} ---> hsp.bits: {} {} ? {}".format(json_data[model_id]["model_name"],hsp.bits,type(hsp.bits), type(true_pass_evalue)))

					blastResults = self.results(blastResults, blast_record.query + " | QUERY: " + alignment.hit_def, perfect, strict , loose, self.include_nudge)

			return blastResults
