from app.settings import *
import csv, glob

class BWT(object):
	"""
	Class to align metagenomic reads to CARD and wildCARD reference using bwa or bowtie2 and 
	provide reports (gene, allele report and read level reports).
	"""
	def __init__(self, aligner, include_wildcard, read_one, read_two, threads, output_file, debug, local_database):
		"""Creates BWT object."""
		self.aligner = aligner
		self.read_one = read_one
		self.read_two = read_two
		self.threads = threads
		self.output_file = output_file

		self.local_database = local_database
		self.db = path
		self.data = data_path
		self.include_wildcard = include_wildcard

		if self.local_database:
			self.db = LOCAL_DATABASE
			self.data = LOCAL_DATABASE

		self.reference_genome = os.path.join(self.data, "card_reference.fasta")

		# index dbs
		self.indecies_directory = os.path.join(self.db,"bwt")
		self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bowtie2"))
		self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bwa"))

		if self.include_wildcard == True:
			self.reference_genome = os.path.join(self.data, "card_wildcard_reference.fasta")
			self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference", "{}".format("bowtie2"))
			self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference", "{}".format("bwa"))

		# outputs
		self.working_directory = os.path.join(os.getcwd())
		self.output_sam_file = os.path.join(self.working_directory, "{}.sam".format(self.output_file))
		self.output_bam_file = os.path.join(self.working_directory, "{}.bam".format(self.output_file))
		self.output_bam_sorted_file = os.path.join(self.working_directory, "{}.sorted.bam".format(self.output_file))
		self.sorted_bam_sorted_file_length_100 = os.path.join(self.working_directory, "{}.sorted.length_100.bam".format(self.output_file))
		self.output_tab = os.path.join(self.working_directory, "{}.txt".format(self.output_file))
		self.output_tab_sequences = os.path.join(self.working_directory, "{}.seqs.txt".format(self.output_file))
		self.output_tab_coverage = os.path.join(self.working_directory, "{}.coverage.txt".format(self.output_file))
		self.output_tab_coverage_all_positions =  os.path.join(self.working_directory, "{}.coverage_all_positions.txt".format(self.output_file))
		self.output_tab_coverage_all_positions_summary = os.path.join(self.working_directory, "{}.coverage_all_positions.summary.txt".format(self.output_file))
		self.model_species_data_type  = os.path.join(self.working_directory, "{}.model_species_data_type.txt".format(self.output_file))
		self.allele_mapping_data_json = os.path.join(self.working_directory, "{}.allele_mapping_data.json".format(self.output_file))
		self.allele_mapping_data_tab = os.path.join(self.working_directory, "{}.allele_mapping_data.txt".format(self.output_file))
		self.gene_mapping_data_tab = os.path.join(self.working_directory, "{}.gene_mapping_data.txt".format(self.output_file))
		self.debug = debug

		if self.debug:
			logger.setLevel(10)

	def __repr__(self):
		"""Returns BWT class full object."""
		return "BWT({}".format(self.__dict__)

	def create_index(self):
		"""
		"""
		if self.aligner == "bowtie2":
			if not os.path.exists(self.index_directory_bowtie2):
				os.makedirs(self.index_directory_bowtie2)
				logger.info("created index at {}".format(self.index_directory_bowtie2))

			os.system("bowtie2-build {reference_genome} {index_directory} --threads {threads}".format(
				index_directory=self.index_directory_bowtie2,
				reference_genome=self.reference_genome,
				threads=self.threads
				)
			)
		else:
			if not os.path.exists(self.index_directory_bwa):
				os.makedirs(self.index_directory_bwa)
				logger.info("created index at {}".format(self.index_directory_bwa))

			os.system("bwa index -p {index_directory} {reference_genome}".format(
				index_directory=self.index_directory_bwa,
				reference_genome=self.reference_genome
				)
			)

	def align_bowtie2_unpaired(self):
		"""
		"""
		os.system("bowtie2 --local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --threads {threads} -x {index_directory} -U {unpaired_reads}  -S {output_sam_file}".format(
			threads=self.threads,
			index_directory=self.index_directory_bowtie2,
			unpaired_reads=self.read_one,
			output_sam_file=self.output_sam_file
			)
		)

	def align_bowtie2(self):
		"""
		"""
		os.system("bowtie2 --local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --threads {threads} -x {index_directory} -1 {read_one} -2 {read_two}  -S {output_sam_file}".format(
			threads=self.threads,
			index_directory=self.index_directory_bowtie2,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=self.output_sam_file
		))

	def align_bwa_single_end_mapping(self):
		"""
		"""
		os.system("bwa mem -M -t {threads} {index_directory} {read_one} > {output_sam_file}".format(
			threads=self.threads,
			index_directory=self.index_directory_bwa,
			read_one=self.read_one,
			output_sam_file=self.output_sam_file
			)
		)

	def align_bwa_paired_end_mapping(self):
		"""
		"""
		os.system("bwa mem -t {threads} {index_directory} {read_one} {read_two} > {output_sam_file}".format(
			threads=self.threads,
			index_directory=self.index_directory_bwa,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=self.output_sam_file
			)
		)

	def convert_sam_to_bam(self):
		"""
		"""
		os.system("samtools view --threads {threads} -b  {input_sam_file} > {output_bam_file}".format(
			threads=self.threads,
			output_bam_file=self.output_bam_file,
			input_sam_file=self.output_sam_file
			)
		)

	def sort_bam(self):
		"""
		"""
		os.system("samtools sort --threads {threads} -T /tmp/aln.sorted -o {sorted_bam_file} {unsorted_bam_file}".format(
			threads=self.threads,
			unsorted_bam_file=self.output_bam_file,
			sorted_bam_file=self.output_bam_sorted_file
			)
		)

	def index_bam(self, bam_file):
		"""
		"""
		os.system("samtools index {input_bam}".format(input_bam=bam_file))

	def extract_alignments_with_length(self, length=100):
		"""
		"""
		os.system("bamtools filter -in {input_bam} -out {output_bam}".format(
			input_bam=self.output_bam_sorted_file,
			output_bam=self.sorted_bam_sorted_file_length_100
		))

	def get_aligned(self):
		"""
		Get stats for aligned reads using 'samtools idxstats'
		"""
		os.system("samtools idxstats {input_bam} | awk '$3 != 0' > {output_tab}".format(
			input_bam=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab
			)
		)

	def get_qname_rname_sequence(self):
		"""
		MAPQ (mapping quality - describes the uniqueness of the alignment, 0=non-unique, >10 probably unique) | awk '$5 > 0'
		"""
		os.system("samtools view --threads {threads} {input_bam} | cut -f 1,2,3,4,5,7 | sort -s -n -k 1,1 > {output_tab}".format(
			threads=self.threads, 
			input_bam=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab_sequences
			)
		)

	def get_coverage(self):
		"""
		"""
		os.system("samtools depth {sorted_bam_file} > {output_tab}".format(
			sorted_bam_file=self.sorted_bam_sorted_file_length_100, 
			output_tab=self.output_tab_coverage
			)
		)

	def get_coverage_all_positions(self):
		"""
		Get converage for all positions using 'genomeCoverageBed'
		BAM file _must_ be sorted by position
		"""
		
		cmd = "genomeCoverageBed -ibam {sorted_bam_file}  -g {reference_genome} > {output_tab}".format(
			sorted_bam_file=self.sorted_bam_sorted_file_length_100,
			reference_genome=self.reference_genome, 
			output_tab=self.output_tab_coverage_all_positions
		)
		os.system(cmd)
		os.system("cat {output_tab} | awk '$2 > 0' | cut -f1,3,4,5 > {output_file}".format(
			output_tab=self.output_tab_coverage_all_positions,
			output_file=self.output_tab_coverage_all_positions_summary
			)
		)

	def get_reads_count(self):
		"""
		Parse tab-delimited file for counts to a dictionary
		"""
		sequences = {}
		with open(self.output_tab, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				sequences[row[0]] = {
					"mapped": row[2],
					"unmapped": row[3],
					"all": format(sum(map(int, [row[2], row[3]])))
				}
			return sequences

	def get_model_details(self):
		"""
		Parse card.json to get each model details
		"""
		models = {}
		categories = {}
		model_name = ""
		try:
			logger.info(os.path.join(self.data, "card.json"))
			with open(os.path.join(self.data, "card.json"), 'r') as jfile:
				data = json.load(jfile)
		except Exception as e:
			logger.error("{}".format(e))
			exit()

		for i in data:
			if i.isdigit():
				categories = {}
				model_name = data[i]["model_name"]
				taxon = []

				if "model_sequences" in data[i]:
					for item in data[i]["model_sequences"]["sequence"]:
						taxa = " ".join(data[i]["model_sequences"]["sequence"][item]["NCBI_taxonomy"]["NCBI_taxonomy_name"].split()[:2])
						if taxa not in taxon:
							taxon.append(taxa)

				for c in data[i]["ARO_category"]:
					if data[i]["ARO_category"][c]["category_aro_class_name"] not in categories.keys():
						categories[data[i]["ARO_category"][c]["category_aro_class_name"]] = []
					if data[i]["ARO_category"][c]["category_aro_name"] not in categories[data[i]["ARO_category"][c]["category_aro_class_name"]]:
						categories[data[i]["ARO_category"][c]["category_aro_class_name"]].append(data[i]["ARO_category"][c]["category_aro_name"])

				models[data[i]["model_id"]] = {
					"model_name": data[i]["model_name"],
					"model_type": data[i]["model_type"],
					"categories": categories,
					"taxon": taxon
				}
		return models

	def get_variant_details(self):
		"""
		Parse tab-delimited to a dictionary for all variants
		"""
		os.system("cat {index_file} | cut -f1,2,6,8,9,10 | sort > {output_file}".format(
			index_file=os.path.join(self.data, "index-for-model-sequences.txt"), 
			output_file=self.model_species_data_type
			)
		)
		variants = {}

		#   0        1        2                   3             4         5
		# ['1280', '1031', 'Escherichia coli', 'ncbi_contig', 'Strict', '99.64']
		# ['438', '1031', 'Klebsiella pneumoniae', 'ncbi_contig', 'Strict', '99.64']
		'''

		1031: {
			'1280': {
				data_type: "ncbi_contig",
				percent_identity: "99.64",
				rgi_criteria: "Strict"
			},
			'438': {
				data_type: "ncbi_contig",
				percent_identity: "99.64",
				rgi_criteria: "Strict"
			},
			'Escherichia coli': [ncbi_contig],
			'Klebsiella pneumoniae': [ncbi_contig],
		}

		'''

		with open(self.model_species_data_type, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				# add new model
				if row[1] not in variants.keys():
					variants.update({
						row[1]: {
							row[2]: [row[3]],
							row[0]: {"data_type":row[3], "rgi_criteria":row[4], "percent_identity":row[5]}
						}
					})
				# update existing model
				else:
					#  check if prev_id is present
					if row[0] not in variants[row[1]].keys():
						# new prev_id
						variants[row[1]].update({
							row[0]: {"data_type":row[3], "rgi_criteria":row[4], "percent_identity":row[5]}
						})
					#  check if pathogen_name is present
					if row[2] not in variants[row[1]].keys():
						# new pathogen_name
						variants[row[1]].update({
							row[2]: [row[3]]
						})
					else:
						if row[3] not in variants[row[1]][row[2]]:
							# add new data type
							variants[row[1]][row[2]].append(row[3])

		return variants

	def get_alignments(self, hit_id, ref_len=0):
		"""
		Parse tab-delimited file into dictionary for mapped reads
		"""
		sequences = []
		with open(self.output_tab_sequences, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				if hit_id == row[2]:
					sequences.append({
						"qname": str(row[0]),
						"flag": str(row[1]),
						"rname": str(row[2]),
						"pos": str(row[3]),
						"mapq": str(row[4]),
						"mrnm": str(row[5])
						})
			return sequences

	def get_coverage_details(self, hit_id):
		"""
		Parse tab-delimited file
		"""
		sequences = {}
		sequences.update({
			hit_id: {
				"covered": 0,
				"uncovered": 0,
				"length": 0
			}
		})

		with open(self.output_tab_coverage_all_positions_summary, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				if hit_id == row[0]:
					sequences[hit_id]["covered"] =  sequences[hit_id]["covered"] + int(row[1])
					sequences[hit_id]["length"] =  int(row[2])
		sequences[hit_id]["uncovered"] =  sequences[hit_id]["length"] - sequences[hit_id]["covered"]

		return sequences

	def get_summary(self):
		"""
		This function uses the following TAB-delimited files:
			<filename>.coverage_all_positions.summary.txt,
			<filename>.txt,
			<filename>.seqs.txt

		------------------------------------------------------------------
		<filename>.txt | samtools idxstats
		------------------------------------------------------------------
		columns:  
				1. reference sequence name
				2. sequence length
				3. # mapped reads
				4. # unmapped reads

		------------------------------------------------------------------
		<filename>.coverage_all_positions.summary.txt | genomeCoverageBed -ibam
		------------------------------------------------------------------
		columns:
				1. chromosome (or entire genome)
				2. depth of coverage from features in input file
				3. number of bases on chromosome (or genome) with depth equal to column 2.
				4. size of chromosome (or entire genome) in base pairs
				5. fraction of bases on chromosome (or entire genome) with depth equal to column 2.

		used 1,3,4,5

		------------------------------------------------------------------
		<filename>.seqs.txt | samtools view
		------------------------------------------------------------------
		columns:
				1.	QNAME	Query template/pair NAME
				2.	FLAG	bitwise FLAG
				3.	RNAME	Reference sequence NAME
				4.	POS		1-based leftmost POSition/coordinate of clipped sequence
				5.	MAPQ	MAPping Quality (Phred-scaled)
				6.	CIGAR	extended CIGAR string
				7.	MRNM	Mate Reference sequence NaMe (`=' if same as RNAME)
				8.	MPOS	1-based Mate POSistion
				9.	TLEN	inferred Template LENgth (insert size)
				10.	SEQ		query SEQuence on the same strand as the reference
				11.	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
				12+. OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

		used 1,2,3,4,5, and 7

		"""
		summary = []
		variants = {}
		models = {}

		logger.info("get_reads_count ...")
		reads = self.get_reads_count()
		logger.info("get_model_details ...")
		models = self.get_model_details()

		if self.include_wildcard:
			logger.info("get_variant_details ...")
			variants = self.get_variant_details()

		mapq_average = 0

		for alignment_hit in reads.keys():
			logger.info(alignment_hit)
			coverage = self.get_coverage_details(alignment_hit)
			model_id = alignment_hit.split("|")[1].split(":")[1]
			cvterm_name = models[model_id]["model_name"]
			model_type = models[model_id]["model_type"]
			resistomes = models[model_id]["categories"]
			alignments = self.get_alignments(alignment_hit)
			mapq_l = []
			mate_pair = []
			for a in alignments:
				mapq_l.append(int(a["mapq"]))
				if a["mrnm"] != "=" and a["mrnm"] not in mate_pair:
					mate_pair.append(a["mrnm"])
					# if cvterm_name.replace(" ", "_") in alignment_hit:
					# 	logger.info("=> {},{}".format(cvterm_name.replace(" ", "_"), alignment_hit))

			if len(mapq_l) > 0:
				mapq_average = sum(mapq_l)/len(mapq_l)

			observed_in_genomes = "no data"
			observed_in_plasmids = "no data"
			prevalence_sequence_id = ""
			observed_data_types = []
			# Genus and species level only (only get first two words)
			observed_in_pathogens = []
			database = "CARD"
			reference_allele_source = "CARD curation"

			if "Prevalence_Sequence_ID" in alignment_hit:
				database = "Resistomes & Variants"
				prevalence_sequence_id = alignment_hit.split("|")[0].split(":")[1]

			if variants:
				if model_id in variants.keys():
					for s in variants[model_id]:
						if s.isdigit() == False:
							observed_in_genomes = "NO"
							observed_in_plasmids = "NO"
							for d in variants[model_id][s]:
								if d not in observed_data_types:
									observed_data_types.append(d)
							if s not in observed_in_pathogens:
								observed_in_pathogens.append(s.replace('"', ""))

					if database != "CARD":
						if "ncbi_chromosome" in observed_data_types:
							observed_in_genomes = "YES"
						if "ncbi_plasmid" in observed_data_types:
							observed_in_plasmids = "YES"

						try:
							reference_allele_source = "In silico {rgi_criteria} {percent_identity}% identity".format(
								rgi_criteria=variants[model_id][prevalence_sequence_id]["rgi_criteria"],
								percent_identity=variants[model_id][prevalence_sequence_id]["percent_identity"],
							)
						except Exception as e:
							reference_allele_source = ""
							logger.warning("missing key with prev_id {} , {}".format(prevalence_sequence_id, e))

				else:
					# provide info from model
					# logger.warning("model not in prev: {}".format(model_id))
					observed_in_pathogens = models[model_id]["taxon"]
			
			# check all clases categories
			if "AMR Gene Family" not in resistomes.keys():
				resistomes["AMR Gene Family"] = []
			if "Drug Class" not in resistomes.keys():
				resistomes["Drug Class"] = []
			if "Resistance Mechanism" not in resistomes.keys():
				resistomes["Resistance Mechanism"] = []

			summary.append({
				"id": alignment_hit,
				"cvterm_name": cvterm_name,
				"model_type": model_type,
				"database": database,
				"reference_allele_source": reference_allele_source,
				"observed_in_genomes": observed_in_genomes,
				"observed_in_plasmids": observed_in_plasmids,
				"observed_in_pathogens": observed_in_pathogens,
				"reads": reads[alignment_hit],
				"alignments": alignments,
				"mapq_average": format(mapq_average,'.2f'),
				"mate_pair": mate_pair,

				"percent_coverage": {
					"covered": format(float(coverage[alignment_hit]["covered"] / coverage[alignment_hit]["length"])*100,'.2f' ),
					"uncovered": format(float(coverage[alignment_hit]["uncovered"] / coverage[alignment_hit]["length"])*100,'.2f')
				},
				"length_coverage": {
					"covered": "{}".format(coverage[alignment_hit]["covered"]),
					"uncovered": "{}".format(coverage[alignment_hit]["uncovered"])
				},
				"reference": {
					"sequence_length": "{}".format(coverage[alignment_hit]["length"])
				},

				"mutation": "N/A",
				"resistomes": resistomes
				,"predicted_pathogen": "N/A"
				})
				# 
			# print(summary)
			# logger.info(json.dumps(summary, indent=2))
			# exit()
		# write json
		with open(self.allele_mapping_data_json, "w") as af:
			af.write(json.dumps(summary,sort_keys=True))

		# wrtie tab-delimited allele_mapping_data
		with open(self.allele_mapping_data_tab, "w") as tab_out:
			writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
			writer.writerow([
							"Reference Sequence",
							"ARO Term",
							"Reference Model Type",
							"Reference DB",
							"Reference Allele Source",
							"Resistomes & Variants: Observed in Genome(s)",
							"Resistomes & Variants: Observed in Plasmid(s)",
							"Observed in Pathogen (CARD & Resistomes)",
							"Completely Mapped Reads",
							"Mapped Reads with Flanking Sequence",
							"All Mapped Reads",
							"Percent Coverage",
							"Length Coverage (bp)",  
							"Average MAPQ (Completely Mapped Reads)",
							"Mate Pair Linkage",
							"Reference Length",
							# "Mutation",
							"AMR Gene Family",
							"Drug Class",
							"Resistance Mechanism"
							# ,"Predicted Pathogen"
							])
			for r in summary:
				writer.writerow([
					r["id"], 
					r["cvterm_name"],
					r["model_type"],
					r["database"],
					r["reference_allele_source"],
					r["observed_in_genomes"],
					r["observed_in_plasmids"],
					"; ".join(r["observed_in_pathogens"]),
					r["reads"]["mapped"], 
					r["reads"]["unmapped"],
					r["reads"]["all"],
					r["percent_coverage"]["covered"],
					r["length_coverage"]["covered"],
					r["mapq_average"],
					"; ".join(r["mate_pair"]),
					r["reference"]["sequence_length"],
					# r["mutation"],
					"; ".join(r["resistomes"]["AMR Gene Family"]),
					"; ".join(r["resistomes"]["Drug Class"]),
					"; ".join(r["resistomes"]["Resistance Mechanism"])
					# ,r["predicted_pathogen"]
				])

		# wrtie tab-delimited gene_mapping_data
		mapping_summary = {}
		alleles_mapped = []
		for r in summary:
			alleles_mapped.append(r["cvterm_name"])
			if r["cvterm_name"] not in mapping_summary.keys():
				# initialise
				mapping_summary[r["cvterm_name"]] = {
					"model_type": [],
					"database": [],
					"alleles_mapped": [],
					"observed_in_genomes": [],
					"observed_in_plasmids": [],
					"observed_in_pathogens": [],
					"mapped": [],
					"unmapped": [],
					"all": [],
					"percent_coverage": [],
					"length_coverage": [],
					"mapq_average": [],
					"mate_pair": [],
					"AMR Gene Family": [],
					"Drug Class": [],
					"Resistance Mechanism": []
				}

				mapping_summary[r["cvterm_name"]]["model_type"].append(r["model_type"])
				mapping_summary[r["cvterm_name"]]["database"].append(r["database"])
				mapping_summary[r["cvterm_name"]]["observed_in_genomes"].append(r["observed_in_genomes"])
				mapping_summary[r["cvterm_name"]]["observed_in_plasmids"].append(r["observed_in_plasmids"])	

				for p in r["observed_in_pathogens"]:
					mapping_summary[r["cvterm_name"]]["observed_in_pathogens"].append(p)

				mapping_summary[r["cvterm_name"]]["mapped"].append(r["reads"]["mapped"])
				mapping_summary[r["cvterm_name"]]["unmapped"].append(r["reads"]["unmapped"])
				mapping_summary[r["cvterm_name"]]["all"].append(r["reads"]["all"])

				mapping_summary[r["cvterm_name"]]["percent_coverage"].append(r["percent_coverage"]["covered"])
				mapping_summary[r["cvterm_name"]]["length_coverage"].append(r["length_coverage"]["covered"])
				mapping_summary[r["cvterm_name"]]["mapq_average"].append(r["mapq_average"])

				for m in r["mate_pair"]:
					# Prevalence_Sequence_ID:12760|ID:1786|Name:mexY
					if m not in ["*"]:
						mapping_summary[r["cvterm_name"]]["mate_pair"].append("{}".format(m.split("|")[2].split(":")[1]))

				for a in r["resistomes"]["AMR Gene Family"]:
					mapping_summary[r["cvterm_name"]]["AMR Gene Family"].append(a)

				for d in r["resistomes"]["Drug Class"]:
					mapping_summary[r["cvterm_name"]]["Drug Class"].append(d)

				for c in r["resistomes"]["Resistance Mechanism"]:
					mapping_summary[r["cvterm_name"]]["Resistance Mechanism"].append(c)

			else:
				if r["model_type"] not in mapping_summary[r["cvterm_name"]]["model_type"]:
					mapping_summary[r["cvterm_name"]]["model_type"].append(r["model_type"])
				if r["database"] not in mapping_summary[r["cvterm_name"]]["database"]:
					mapping_summary[r["cvterm_name"]]["database"].append(r["database"])
				if r["observed_in_genomes"] not in mapping_summary[r["cvterm_name"]]["observed_in_genomes"]:
					mapping_summary[r["cvterm_name"]]["observed_in_genomes"].append(r["observed_in_genomes"])
				if r["observed_in_plasmids"] not in mapping_summary[r["cvterm_name"]]["observed_in_plasmids"]:
					mapping_summary[r["cvterm_name"]]["observed_in_plasmids"].append(r["observed_in_plasmids"])	

				# loop thru array
				for p in r["observed_in_pathogens"]:
					if p not in mapping_summary[r["cvterm_name"]]["observed_in_pathogens"]:
						mapping_summary[r["cvterm_name"]]["observed_in_pathogens"].append(p)	

				mapping_summary[r["cvterm_name"]]["mapped"].append(r["reads"]["mapped"])
				mapping_summary[r["cvterm_name"]]["unmapped"].append(r["reads"]["unmapped"])
				mapping_summary[r["cvterm_name"]]["all"].append(r["reads"]["all"])	

				mapping_summary[r["cvterm_name"]]["percent_coverage"].append(r["percent_coverage"]["covered"])
				mapping_summary[r["cvterm_name"]]["length_coverage"].append(r["length_coverage"]["covered"])
				mapping_summary[r["cvterm_name"]]["mapq_average"].append(r["mapq_average"])

				for m in r["mate_pair"]:
					if m not in ["*"]:
						mapping_summary[r["cvterm_name"]]["mate_pair"].append("{}".format(m.split("|")[2].split(":")[1]))

				for a in r["resistomes"]["AMR Gene Family"]:
					if a not in mapping_summary[r["cvterm_name"]]["AMR Gene Family"]:
						mapping_summary[r["cvterm_name"]]["AMR Gene Family"].append(a)

				for d in r["resistomes"]["Drug Class"]:
					if d not in mapping_summary[r["cvterm_name"]]["Drug Class"]:
						mapping_summary[r["cvterm_name"]]["Drug Class"].append(d)

				for c in r["resistomes"]["Resistance Mechanism"]:
					if c not in mapping_summary[r["cvterm_name"]]["Resistance Mechanism"]:
						mapping_summary[r["cvterm_name"]]["Resistance Mechanism"].append(c)


		with open(self.gene_mapping_data_tab, "w") as tab_out:
			writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
			writer.writerow([
							"ARO Term",
							"Reference Model Type",
							"Reference DB",
							"Alleles Mapped",
							"Resistomes & Variants: Observed in Genome(s)",
							"Resistomes & Variants: Observed in Plasmid(s)",
							"Observed Pathogen (CARD & Resistomes)",
							"Completely Mapped Reads",
							"Mapped Reads with Flanking Sequence",
							"All Mapped Reads",
							"Average Percent Coverage",
							"Average Length Coverage (bp)",  
							"Average MAPQ (Completely Mapped Reads)",
							"Mate Pair Linkage (# reads)",
							"AMR Gene Family",
							"Drug Class",
							"Resistance Mechanism"
							])
			am = { item:alleles_mapped.count(item) for item in alleles_mapped }
			for i in mapping_summary:
				observed_in_genomes = "NO"
				observed_in_plasmids = "NO"

				if "YES" in mapping_summary[i]["observed_in_genomes"]:
					observed_in_genomes = "YES"
				elif "no data" in mapping_summary[i]["observed_in_genomes"]:
					observed_in_genomes = "no data"

				if "YES" in mapping_summary[i]["observed_in_plasmids"]:
					observed_in_plasmids = "YES"
				elif "no data" in mapping_summary[i]["observed_in_plasmids"]:
					observed_in_plasmids = "no data"

				average_percent_coverage = 0
				average_length_coverage = 0
				average_mapq  = 0
		
				if len(mapping_summary[i]["percent_coverage"]) > 0:
					average_percent_coverage = sum(map(float,mapping_summary[i]["percent_coverage"]))/len(mapping_summary[i]["percent_coverage"])

				if len(mapping_summary[i]["length_coverage"]) > 0:
					average_length_coverage = sum(map(float,mapping_summary[i]["length_coverage"]))/len(mapping_summary[i]["length_coverage"])

				if len(mapping_summary[i]["mapq_average"]) > 0:
					average_mapq = sum(map(float,mapping_summary[i]["mapq_average"]))/len(mapping_summary[i]["mapq_average"])

				mate_pairs = []
				mp = { item:mapping_summary[i]["mate_pair"].count(item) for item in mapping_summary[i]["mate_pair"]}
				for k in mp:
					if k != i.replace(" ", "_"):
						mate_pairs.append("{} ({})".format(k,mp[k]))

				writer.writerow([
					i,
					"; ".join(mapping_summary[i]["model_type"]),
					"; ".join(mapping_summary[i]["database"]),
					am[i],
					observed_in_genomes,
					observed_in_plasmids,
					"; ".join(mapping_summary[i]["observed_in_pathogens"]),
					format(sum(map(float,mapping_summary[i]["mapped"])),'.2f'),
					format(sum(map(float,mapping_summary[i]["unmapped"])),'.2f'),
					format(sum(map(float,mapping_summary[i]["all"])),'.2f'),
					format(average_percent_coverage,'.2f'),
					format(average_length_coverage,'.2f'),
					format(average_mapq,'.2f'),
					"; ".join(mate_pairs),
					"; ".join(mapping_summary[i]["AMR Gene Family"]),
					"; ".join(mapping_summary[i]["Drug Class"]),
					"; ".join(mapping_summary[i]["Resistance Mechanism"])
				])

	def check_index(self):
		"""
		Check if index exists for a given reference fasta file.
		"""
		# check if we have db
		files = [os.path.basename(x) for x in glob.glob(os.path.join(self.indecies_directory,"*"))]
		logger.info(json.dumps(files, indent=2))
		if self.aligner == "bowtie2":
			if (("bowtie2.1.bt2" in files) and \
				("bowtie2.2.bt2" in files) and \
				("bowtie2.3.bt2" in files) and \
				("bowtie2.4.bt2" in files) and \
				("bowtie2.rev.1.bt2" in files) and \
				("bowtie2.rev.2.bt2" in files)) == False:
				# create index and save results in ./db from reference genome: (.fasta)
				self.create_index()
			else:
				logger.info("{} index already exists".format(self.aligner))
		else:
			if (("bwa.amb" in  files) and \
				("bwa.ann" in  files) and \
				("bwa.bwt" in  files) and \
				("bwa.pac" in  files) and \
				("bwa.sa" in files)) == False:
				# create index and save results in ./db from reference genome: (.fasta)
				self.create_index()
			else:
				logger.info("{} index already exists".format(self.aligner))


	def run(self):
		"""
		Align reads to reference genomes and report
		"""
		# print args
		logger.info(json.dumps(self.__dict__, indent=2))

		# index database
		logger.info("index database")
		self.check_index()

	    # align
		logger.info("align using {}".format(self.aligner))
		if self.aligner == "bowtie2":
			if self.read_two == None:
				self.align_bowtie2_unpaired()
			else:
				self.align_bowtie2()
		else:
			if self.read_two == None:
				self.align_bwa_single_end_mapping()
			else:
				self.align_bwa_paired_end_mapping()

		# convert SAM file to BAM file
		logger.info("convert SAM file to BAM file")
		self.convert_sam_to_bam()

		# sort BAM file
		logger.info("sort BAM file")
		self.sort_bam()

		# index BAM file
		logger.info("index BAM file")
		self.index_bam(bam_file=self.output_bam_sorted_file)

		# only extract alignment of specific length 
		logger.info("only extract alignment of specific length")
		self.extract_alignments_with_length()

		# index filtered BAM file
		logger.info("index filtered BAM file")
		self.index_bam(bam_file=self.sorted_bam_sorted_file_length_100)

		# pull alligned
		logger.info("pull alligned")
		self.get_aligned()

		# pull qname, rname and sequence
		logger.info("pull qname, rname and sequence")
		self.get_qname_rname_sequence()

		# get coverage
		logger.info("get coverage")
		self.get_coverage()

		# get coverage for all positions
		logger.info("get coverage for all positions")
		self.get_coverage_all_positions()

		# get summary
		logger.info("get summary")
		self.get_summary()

		logger.info("Done.")













