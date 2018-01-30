from app.settings import *

class Database(object):
	"""Class to create BLAST databases from a card.json file."""
	def __init__(self):
		"""Creates Database object."""

	def __repr__(self):
		"""Returns Database class full object."""
		return "Database({}".format(self.__dict__)

	def build_databases(self):
		"""Build BLAST and DIAMOND databases."""
		self.write_fasta_from_json()
		self.make_blast_database()
		self.make_diamond_database()
		self.write_fasta_from_json_rna()

	def make_blast_database(self):
		"""Build BLAST database from a FASTA file."""
		if os.path.isfile(path+"proteindb.fsa") == True and os.path.exists(path+"proteindb.fsa") == True  \
		   and os.path.exists(path+"protein.db.phr") == True and os.path.exists(path+"protein.db.pin") == True \
		   and os.path.exists(path+"protein.db.psq") == True:
		   logger.info("blast DB exists")
		   pass
		else:
			logger.info("create blast DB.")
			os.system('makeblastdb -in '+path+'proteindb.fsa -dbtype prot -out '+path+'protein.db 2>&1 >> ' + file_handler.baseFilename)

	def make_diamond_database(self):
		"""Build DIAMOND database from a FASTA file."""
		if os.path.isfile(path+"proteindb.fsa") == True and os.path.exists(path+"proteindb.fsa") == True \
			and os.path.exists(path+"protein.db.dmnd") == True:
			logger.info("diamond DB exists")
			pass
		else:
			logger.info("create diamond DB.")
			os.system('diamond makedb --quiet --in '+path+'proteindb.fsa --db '+path+'protein.db 2>&1 >> ' + file_handler.baseFilename)

	def make_custom_db(self, in_file, out_file, db_type="nucl", program="blast"):
		if program == 'blast':
			os.system('makeblastdb -in {in_file}  -dbtype {db_type} -out {out_file} 2>&1 >> \
				{log_file}'.format( in_file=in_file, db_type=db_type, out_file=out_file, log_file=file_handler.baseFilename))
		else:
			exit("Only NCBI BLAST is supported.")

	def write_fasta_from_json(self):
		"""Creates a fasta file from card.json file."""
		if os.path.isfile('%sproteindb.fsa' % (path)):
			# logger.info("Database already exists.")
			return
		else:
			with open('%scard.json' % (data_path), 'r') as jfile:
				j = json.load(jfile)

			with open('%sproteindb.fsa' % (path), 'w') as fout:
				for i in j:
					if i.isdigit():
		            	# model_type: protein homolog model
						if j[i]['model_type_id'] == '40292':
							try:
								pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
							except KeyError:
								logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
							else:
								for seq in j[i]['model_sequences']['sequence']:
									fout.write('>%s_%s | model_type_id: 40292 | pass_bitscore: %s | %s\n' % (i, seq, pass_bit_score, j[i]['ARO_name']))
									fout.write('%s\n' %(j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))

		            	# model_type: protein variant model
						elif j[i]["model_type_id"] == "40293":
							try:
								pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
							except KeyError:
								logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
							else:
								snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
								for seq in j[i]['model_sequences']['sequence']:
									fout.write('>%s_%s | model_type_id: 40293 | pass_bit_score: %s | SNP: %s | %s\n' \
										% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
									fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))

		            	# model_type: protein overexpression model
						elif j[i]["model_type_id"] == "41091":
							try:
								pass_bit_score = j[i]["model_param"]["blastp_bit_score"]["param_value"]
							except KeyError:
								logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
							else:
								snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
								for seq in j[i]['model_sequences']['sequence']:
									fout.write('>%s_%s | model_type_id: 41091 | pass_bit_score: %s | SNP: %s | %s\n' \
										% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
									fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))


	def write_fasta_from_json_rna(self):
		snpList_16s = []
		snpList_23s = []
		"""Creates a fasta file for 16S and 23S data from card.json file."""
		if os.path.isfile('%srnadb.fsa' % (path)):
			# logger.info("RNA database already exists.")
			return
		else:
			with open('%scard.json' % (data_path), 'r') as jfile:
				j = json.load(jfile)

			with open('%srnadb.fsa' % (path), 'w') as fout:
				for i in j:
					if i.isdigit():
						# model_type: ribosomal RNA model
						if j[i]["model_type_id"] == "40295":
							try:
								pass_bit_score = j[i]['model_param']['blastn_bit_score']['param_value']
							except KeyError:
								logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
									% (j[i]['model_id'], j[i]['model_name']))
								logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
							else:
								snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]

								for s in snpList:
									if "16S" in j[i]['ARO_name']:
										if s not in snpList_16s:
											snpList_16s.append(s)
									if "23S" in j[i]['ARO_name']:
										if s not in snpList_23s:
											snpList_23s.append(s)

								for seq in j[i]['model_sequences']['sequence']:
								
									
									if j[i]['model_sequences']['sequence'][seq]['dna_sequence']['strand'] == "-":
										basecomplement = self.complementary_strand(j[i]['model_sequences']['sequence'][seq]['dna_sequence']['sequence'])
							
										fout.write('>%s_%s | model_type_id: 40295 | pass_bit_score: %s | SNP: %s | %s\n' \
										% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
										fout.write('%s\n' % (basecomplement))

									else:
									
										fout.write('>%s_%s | model_type_id: 40295 | pass_bit_score: %s | SNP: %s | %s\n' \
										% (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
										fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['dna_sequence']['sequence']))

		# write snps to file
		with open('%s16s_rRNA.txt' % (path), 'w') as f16s:
			snpList_16s.sort()
			f16s.write(',\n'.join(snpList_16s))

		with open('%s23s_rRNA.txt' % (path), 'w') as f23s:
			snpList_23s.sort()
			f23s.write(',\n'.join(snpList_23s))


	def complementary_strand(self, strand):
		'''Takes a DNA strand string and returns its opposite base pair match.'''
		self.trans = { "T": "A", "A": "T", "G": "C", "C": "G" , "N": "N", "M":"K", "K":"M", \
					   "R":"Y", "Y":"R", "S":"S", "W":"W", "B":"V", "V":"B", "H":"D", "D":"H"}
		complement = []

		for base in strand: 
			complement.append(self.trans[base])
		
		complement_seq = ''.join(complement)
		return complement_seq



