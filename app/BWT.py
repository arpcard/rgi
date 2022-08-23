from app.settings import *
import csv, glob
from multiprocessing import Pool
import time
import statistics
from Bio import SeqIO, Seq
import gzip
# pip3 install "dask[dataframe]"
import dask.dataframe as dd
import subprocess

class BWT(object):
	"""
	Class to align metagenomic reads to CARD and wildCARD reference using bwa or bowtie2 and
	provide reports (gene, allele report and read level reports).
	"""
	def __init__(self, aligner, include_wildcard, include_baits, read_one, read_two, threads, output_file, debug, clean, local_database, mapq, mapped, coverage, include_other_models):
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
		self.include_baits = include_baits

		self.mapq = mapq
		self.mapped = mapped
		self.coverage = coverage
		self.include_other_models = include_other_models

		self.debug = debug
		if self.debug:
			logger.setLevel(10)

		if self.local_database:
			self.db = LOCAL_DATABASE
			self.data = LOCAL_DATABASE

		# index dbs
		self.indecies_directory = os.path.join(self.db,"bwt")

		"""
		Files:

		card_reference.fasta == CARD[homolog]
		card_reference_all.fasta == CARD[homolog,rRNA,overexpression,variants,knockout]
		card_wildcard_reference.fasta == CARD[homolog] + WILD[homolog]
		card_wildcard_reference_all.fasta == CARD[homolog,rRNA,overexpression,variants,knockout] + WILD[homolog,rRNA,overexpression,variants,knockout]
		"""
		# card canonical
		# CARD[homolog]
		# logger.info("card canonical")
		# self.reference_genome = os.path.join(self.data, "card_reference.fasta")
		# self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bowtie2"))
		# self.index_directory_kma = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("kma"))
		# self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bwa"))

		if self.include_baits == True:
			logger.info("baits")
			self.reference_genome_baits = os.path.join(self.data, "baits_reference.fasta")
			self.index_directory_bowtie2_baits = os.path.join(self.db, self.indecies_directory, "baits_reference", "{}".format("bowtie2_baits"))
			self.index_directory_bwa_baits = os.path.join(self.db, self.indecies_directory, "baits_reference", "{}".format("bwa_baits"))

		# card canonical and variants
		if self.include_wildcard == True and self.include_other_models == False:
			logger.info("CARD[homolog] + WILD[homolog]")
			self.reference_genome = os.path.join(self.data, "card_wildcard_reference.fasta")
			self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference", "{}".format("bowtie2"))
			self.index_directory_kma = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference", "{}".format("kma"))
			self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference", "{}".format("bwa"))
		elif self.include_wildcard == False and self.include_other_models == True:
			logger.info("CARD[homolog,rRNA,overexpression,variants,knockout]")
			self.reference_genome = os.path.join(self.data, "card_reference_all.fasta")
			self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_reference_all", "{}".format("bowtie2"))
			self.index_directory_kma = os.path.join(self.db, self.indecies_directory, "card_reference_all", "{}".format("kma"))
			self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_reference_all", "{}".format("bwa"))
		elif self.include_wildcard == True and self.include_other_models == True:
			logger.info("CARD[homolog,rRNA,overexpression,variants,knockout] + WILD[homolog,rRNA,overexpression,variants,knockout]")
			self.reference_genome = os.path.join(self.data, "card_wildcard_reference_all.fasta")
			self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference_all", "{}".format("bowtie2"))
			self.index_directory_kma = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference_all", "{}".format("kma"))
			self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_wildcard_reference_all", "{}".format("bwa"))
		else:
			logger.info("CARD[homolog]")
			self.reference_genome = os.path.join(self.data, "card_reference.fasta")
			self.index_directory_bowtie2 = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bowtie2"))
			self.index_directory_kma = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("kma"))
			self.index_directory_bwa = os.path.join(self.db, self.indecies_directory, "card_reference", "{}".format("bwa"))

		# outputs
		self.working_directory = os.path.join(os.getcwd())
		self.output_sam_file = os.path.join(self.working_directory, "{}.temp.sam".format(self.output_file))
		self.output_sam_file_baits = os.path.join(self.working_directory, "{}.baits.temp.sam".format(self.output_file))
		self.output_bam_file = os.path.join(self.working_directory, "{}.temp.bam".format(self.output_file))
		self.output_bam_file_baits = os.path.join(self.working_directory, "{}.baits.temp.bam".format(self.output_file))
		self.output_bam_sorted_file = os.path.join(self.working_directory, "{}.sorted.temp.bam".format(self.output_file))
		self.sorted_bam_sorted_file_length_100 = os.path.join(self.working_directory, "{}.sorted.length_100.bam".format(self.output_file))

		self.unmapped = os.path.join(self.working_directory, "{}.unmapped.temp.bam".format(self.output_file))
		self.mapped = os.path.join(self.working_directory, "{}.mapped.temp.bam".format(self.output_file))
		self.mapping_overall_stats = os.path.join(self.working_directory, "{}.overall_mapping_stats.txt".format(self.output_file))
		self.mapping_artifacts_stats = os.path.join(self.working_directory, "{}.artifacts_mapping_stats.txt".format(self.output_file))
		self.mapping_reference_stats = os.path.join(self.working_directory, "{}.reference_mapping_stats.txt".format(self.output_file))
		self.mapping_baits_stats = os.path.join(self.working_directory, "{}.baits_mapping_stats.txt".format(self.output_file))

		self.baits_reads_count = os.path.join(self.working_directory, "{}.baits_reads_count.temp.txt".format(self.output_file))
		self.reads_baits_count = os.path.join(self.working_directory, "{}.reads_baits_count.temp.txt".format(self.output_file))
		self.aro_term_reads = os.path.join(self.working_directory, "{}.aro_term_reads.temp.txt".format(self.output_file))

		self.output_tab = os.path.join(self.working_directory, "{}.temp.txt".format(self.output_file))
		self.output_tab_sequences = os.path.join(self.working_directory, "{}.seqs.temp.txt".format(self.output_file))
		self.output_tab_coverage = os.path.join(self.working_directory, "{}.coverage.temp.txt".format(self.output_file))
		self.output_tab_coverage_all_positions =  os.path.join(self.working_directory, "{}.coverage_all_positions.temp.txt".format(self.output_file))
		self.output_tab_coverage_all_positions_summary = os.path.join(self.working_directory, "{}.coverage_all_positions.summary.temp.txt".format(self.output_file))
		self.model_species_data_type  = os.path.join(self.working_directory, "{}.model_species_data_type.temp.txt".format(self.output_file))
		self.allele_mapping_data_json = os.path.join(self.working_directory, "{}.allele_mapping_data.json".format(self.output_file))
		self.allele_mapping_data_tab = os.path.join(self.working_directory, "{}.allele_mapping_data.txt".format(self.output_file))
		self.gene_mapping_data_tab = os.path.join(self.working_directory, "{}.gene_mapping_data.txt".format(self.output_file))

		self.baits_mapping_data_tab = os.path.join(self.working_directory, "{}.baits_mapping_data.temp.txt".format(self.output_file))
		self.baits_mapping_data_json = os.path.join(self.working_directory, "{}.baits_mapping_data.temp.json".format(self.output_file))
		self.reads_mapping_data_json = os.path.join(self.working_directory, "{}.reads_mapping_data.temp.json".format(self.output_file))

		# map baits to complete genes
		self.baits_card_sam = os.path.join(self.working_directory, "{}.baits_card.temp.sam".format(self.output_file))
		self.baits_card_bam = os.path.join(self.working_directory, "{}.baits_card.temp.bam".format(self.output_file))
		self.baits_card_tab = os.path.join(self.working_directory, "{}.baits_card.temp.txt".format(self.output_file))
		self.baits_card_json = os.path.join(self.working_directory, "{}.baits_card.temp.json".format(self.output_file))
		self.baits_card_data_tab = os.path.join(self.working_directory, "{}.baits_card_data.temp.txt".format(self.output_file))
		self.card_baits_reads_count_json = os.path.join(self.working_directory, "{}.card_baits_reads_count.temp.json".format(self.output_file))


		self.clean = clean


	def __repr__(self):
		"""Returns BWT class full object."""
		return "BWT({}".format(self.__dict__)


	def clean_files(self):
		"""Cleans temporary files."""
		if self.clean == True:
			basename_output_file = os.path.splitext(os.path.basename(self.output_file))[0]
			logger.info("Cleaning up temporary files...{}".format(basename_output_file))
			# clean working_directory
			self.clean_directory(self.working_directory, basename_output_file)
			d_name, f_name = os.path.split(self.output_file)
			# clean destination_directory
			self.clean_directory(d_name, basename_output_file)
		else:
			logger.info("Clean up skipped.")

	def clean_directory(self, directory, basename_output_file):
		"""Cleans files in directory."""
		logger.info(directory)
		files = glob.glob(os.path.join(directory, "*"))
		for f in files:
			if os.path.basename(self.output_file) in f and ".temp" in f and os.path.isfile(f):
				self.remove_file(f)

	def remove_file(self, f):
		"""Removes file."""
		if os.path.exists(f):
			try:
				logger.info("Removed file: {}".format(f))
				os.remove(f)
			except Exception as e:
				raise e
		else:
			logger.warning("Missing file: {}".format(f))


	def create_index(self, index_directory, reference_genome):
		"""
		Create kma, bwa or bowtie2 index for reference genome
		"""
		if self.aligner == "bowtie2":
			if not os.path.exists(index_directory):
				os.makedirs(index_directory)
				logger.info("created index at {}".format(index_directory))

			os.system("bowtie2-build --quiet {reference_genome} {index_directory} --threads {threads}".format(
				index_directory=index_directory,
				reference_genome=reference_genome,
				threads=self.threads
				)
			)
		elif self.aligner == "kma":
			if not os.path.exists(index_directory):
				os.makedirs(index_directory)
				logger.info("created index at {}".format(index_directory))
			os.system("kma index -i {reference_genome} -o {index_directory}".format(
				index_directory=index_directory,
				reference_genome=reference_genome
				)
			)
		else:
			if not os.path.exists(index_directory):
				os.makedirs(index_directory)
				logger.info("created index at {}".format(index_directory))

			os.system("bwa index -p {index_directory} {reference_genome}".format(
				index_directory=index_directory,
				reference_genome=reference_genome
				)
			)

	def align_kma(self, reference_genome, index_directory, output_sam_file):
		"""
		Align paired reads to reference database selected (i.e card, wildcard etc) using kma
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)

		logger.info("align reads -ipe {} {} to {}".format(self.read_one, self.read_two, reference_genome))

		cmd = "kma -mem_mode -ex_mode -1t1 -vcf -ipe {read_one} {read_two} -t {threads} -t_db {index_directory} -o {output_sam_file}.temp -sam > {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=output_sam_file
		)
		os.system(cmd)

	def align_kma_interleaved(self, reference_genome, index_directory, output_sam_file):
		"""
		Align interleaved reads to reference database selected (i.e card, wildcard etc) using kma
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)

		logger.info("align reads -ipe {} {} to {}".format(self.read_one, self.read_two, reference_genome))

		cmd = "kma -mem_mode -ex_mode -1t1 -vcf -int {read_one} -t {threads} -t_db {index_directory} -o {output_sam_file}.temp -sam > {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=output_sam_file
		)
		os.system(cmd)

	def align_bowtie2_unpaired(self, reference_genome, index_directory, output_sam_file):
		"""
		Align unpaired reads to reference database selected (i.e card, wildcard etc) using bowtie2
		"""

		self.check_index(index_directory=index_directory, reference_genome=reference_genome)

		cmd = "bowtie2 --very-sensitive-local --threads {threads} -x {index_directory} -U {unpaired_reads}  -S {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			unpaired_reads=self.read_one,
			output_sam_file=output_sam_file
		)

		os.system(cmd)

	def align_bowtie2(self, reference_genome, index_directory, output_sam_file):
		"""
		Align paired reads to reference database selected (i.e card, wildcard etc) using bowtie2
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)

		logger.info("align reads -1 {} -2 {} to {}".format(self.read_one, self.read_two, reference_genome))

		cmd = "bowtie2 --quiet --very-sensitive-local --threads {threads} -x {index_directory} -1 {read_one} -2 {read_two}  -S {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=output_sam_file
		)

		os.system(cmd)

	def align_bowtie2_baits_to_genes(self, reference_genome, index_directory, output_sam_file):
		"""
		Align baits to genes
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)

		logger.info("align baits -f {} to complete genes in {}".format(self.reference_genome_baits, reference_genome))

		cmd = "bowtie2 --quiet --very-sensitive-local --threads {threads} -x {index_directory} -f {unpaired_reads}   -S {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			unpaired_reads=self.reference_genome_baits,
			output_sam_file=output_sam_file
		)

		os.system(cmd)

	def align_bwa_single_end_mapping(self, reference_genome, index_directory, output_sam_file):
		"""
		Align unpaired reads to reference database selected (i.e card, wildcard etc) using bwa
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)
		os.system("bwa mem -M -t {threads} {index_directory} {read_one} > {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			read_one=self.read_one,
			output_sam_file=output_sam_file
			)
		)

	def align_bwa_paired_end_mapping(self, reference_genome, index_directory, output_sam_file):
		"""
		Align paired reads to reference database selected (i.e card, wildcard etc) using bwa
		"""
		self.check_index(index_directory=index_directory, reference_genome=reference_genome)
		os.system("bwa mem -t {threads} {index_directory} {read_one} {read_two} > {output_sam_file}".format(
			threads=self.threads,
			index_directory=index_directory,
			read_one=self.read_one,
			read_two=self.read_two,
			output_sam_file=output_sam_file
			)
		)

	def convert_sam_to_bam(self, input_sam_file, output_bam_file):
		"""
		Convert sam file to bam file
		"""
		os.system("samtools view --threads {threads} -b  {input_sam_file} > {output_bam_file}".format(
			threads=self.threads,
			output_bam_file=output_bam_file,
			input_sam_file=input_sam_file
			)
		)

	def sort_bam(self):
		"""
		Sort bam file
		"""
		os.system("samtools sort --threads {threads} -T {output_file}.sorted -o {sorted_bam_file} {unsorted_bam_file}".format(
			threads=self.threads,
			output_file=self.output_file,
			unsorted_bam_file=self.output_bam_file,
			sorted_bam_file=self.output_bam_sorted_file
			)
		)

	def index_bam(self, bam_file):
		"""
		Index input bam file using 'samtools index' function
		"""
		os.system("samtools index {input_bam}".format(input_bam=bam_file))

	def extract_alignments_with_length(self, length=10, map_quality=2):
		"""
		Get alignments from bam file using length as a filter (default length=10, map_quality=2)
		TODO:: add filters (mapped, length, coverage and map_quality)
		"""
		cmd="bamtools filter -in {input_bam} -out {output_bam}".format(
			input_bam=self.output_bam_sorted_file,
			output_bam=self.sorted_bam_sorted_file_length_100,
			length=length,
			map_quality=map_quality
		)
		# logger.info(cmd)
		subprocess.run(["bamtools","filter","-in",self.output_bam_sorted_file,"-out",self.sorted_bam_sorted_file_length_100])

	def get_aligned(self):
		"""
		Get stats for aligned reads using 'samtools idxstats' function
		"""
		cmd = "samtools idxstats {input_bam} > {output_tab}".format(
			input_bam=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab
			)
		# logger.debug(cmd)
		input_bam=self.sorted_bam_sorted_file_length_100
		output_tab=self.output_tab
		subprocess.run(f"samtools idxstats {input_bam} > {output_tab}", shell=True, check=True)

	def get_qname_rname_sequence(self):
		"""
		MAPQ (mapping quality - describes the uniqueness of the alignment, 0=non-unique, >10 probably unique) | awk '$5 > 0'
		"""
		cmd="samtools view --threads {threads} {input_bam} | cut -f 1,2,3,4,5,7 | sort -s -n -k 1,1 > {output_tab}".format(
			threads=self.threads,
			input_bam=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab_sequences
			)
		# logger.debug(cmd)
		# os.system(cmd)
		threads=self.threads
		input_bam=self.sorted_bam_sorted_file_length_100
		output_tab=self.output_tab_sequences
		subprocess.run(f"samtools view --threads {threads} {input_bam} | cut -f 1,2,3,4,5,7 | sort -s -n -k 1,1 > {output_tab}", shell=True, check=True)

	def get_coverage(self):
		"""
		Get coverage using 'samtools depth' function and write outputs to a tab-delimited file
		"""
		cmd="samtools depth {sorted_bam_file} > {output_tab}".format(
			sorted_bam_file=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab_coverage
			)
		# logger.debug(cmd)
		# os.system(cmd)
		sorted_bam_file=self.sorted_bam_sorted_file_length_100
		output_tab=self.output_tab_coverage
		subprocess.run(f"samtools depth {sorted_bam_file} > {output_tab}", shell=True, check=True)

	def get_coverage_all_positions(self):
		"""
		Get coverage for all positions using 'genomeCoverageBed -ibam' function
		BAM file _must_ be sorted by position
		By default, bedtools genomecov will compute a histogram of coverage for the genome file provided. The default output format is as follows:

		1. chromosome (or entire genome)
		2. depth of coverage from features in input file
		3. number of bases on chromosome (or genome) with depth equal to column 2.
		4. size of chromosome (or entire genome) in base pairs
		5. fraction of bases on chromosome (or entire genome) with depth equal to column 2.

		"""
		# TODO:: This command is slow, it needs improvement, maybe explore -bg and -bga options
		cmd = "bedtools genomecov -ibam {sorted_bam_file} > {output_tab}".format(
			sorted_bam_file=self.sorted_bam_sorted_file_length_100,
			output_tab=self.output_tab_coverage_all_positions
		)
		# logger.debug(cmd)
		# os.system(cmd)
		sorted_bam_file=self.sorted_bam_sorted_file_length_100
		output_tab=self.output_tab_coverage_all_positions
		subprocess.run(f"bedtools genomecov -ibam {sorted_bam_file} > {output_tab}", shell=True, check=True)

		# os.system("cat {output_tab} | awk '$2 > 0' | cut -f1,3,4,5 > {output_file}".format(
		# 	output_tab=self.output_tab_coverage_all_positions,
		# 	output_file=self.output_tab_coverage_all_positions_summary
		# 	)
		# )
		output_tab=self.output_tab_coverage_all_positions
		output_file=self.output_tab_coverage_all_positions_summary
		subprocess.run(f"cat {output_tab} | awk '$2 > 0' | cut -f1,3,4,5 > {output_file}", shell=True, check=True)

	def get_baits_count(self):
		"""
		TODO:: Get baits count
		"""
		pass

	def get_reads_count(self):
		"""
		Parse tab-delimited file for read counts to a dictionary
		"""
		sequences = {}

		with open(self.output_tab, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				if int(row[2]) > 0:
					sequences[row[0]] = {
						"mapped": row[2],
						"unmapped": row[3],
						"all": format(sum(map(int, [row[2], row[3]])))
					}


		# write file reference_stats
		with open(self.mapping_reference_stats, "w") as out:
			out.write("**********************************************\n")
			out.write("Stats for Reference: \n")
			out.write("**********************************************\n")
			out.write("\n")
			out.write("how many reference terms (or ARO terms): {}\n\n".format(len(sequences.keys())))

			return sequences

	def get_model_details(self):
		"""
		Parse card.json to get each model details
		"""
		models = {}
		models_by_accession = {}
		try:
			with open(os.path.join(self.data, "card.json"), 'r') as jfile:
				data = json.load(jfile)
		except Exception as e:
			logger.error("{}".format(e))
			exit()

		for i in data:
			if i.isdigit():
				categories = {}
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
					"model_id": data[i]["model_id"],
					"ARO_accession": data[i]["ARO_accession"],
					"model_name": data[i]["model_name"],
					"model_type": data[i]["model_type"],
					"categories": categories,
					"taxon": taxon
				}

				models_by_accession[data[i]["ARO_accession"]] = {
					"model_id": data[i]["model_id"],
					"ARO_accession": data[i]["ARO_accession"],
					"model_name": data[i]["model_name"],
					"model_type": data[i]["model_type"],
					"categories": categories,
					"taxon": taxon
				}
		return models, models_by_accession

	def get_variant_details(self):
		"""
		Parse tab-delimited to a dictionary for all variants
		"""
		# os.system("cat {index_file} | cut -f1,2,6,7,8,9,10 | sort > {output_file}".format(
		# 	index_file=os.path.join(self.data, "index-for-model-sequences.txt"),
		# 	output_file=self.model_species_data_type
		# 	)
		# )

		df = dd.read_csv(os.path.join(self.data, "index-for-model-sequences.txt"),
			sep="\t",assume_missing=True,usecols=["prevalence_sequence_id","model_id","species_name","ncbi_accession","data_type","rgi_criteria","percent_identity"])
		# , header=None
		df = df.astype({"prevalence_sequence_id": int, "model_id": int})
		df.to_csv(filename=self.model_species_data_type, single_file=True, mode="w", sep="\t", header_first_partition_only=True, index=False)

		variants = {}
		# prevalence_sequence_id	model_id	accession 			species_name			data_type		rgi_criteria	percent_identity
		# 10687						2882		FBHN01				Campylobacter coli		ncbi_contig		Strict			67.71
		#   0        				1        	2                   3             			4         		5       		6
		'''
		1031: { # model_id
			'1280': { # prevalence_sequence_id
				0: { # ncbi accession
					data_type: "ncbi_contig",
					percent_identity: "99.64",
					rgi_criteria: "Strict",
					species_name: 'Escherichia coli'
				}
			},
			'438': { # prevalence_sequence_id
				0: { # ncbi accession
					data_type: "ncbi_contig",
					percent_identity: "99.64",
					rgi_criteria: "Strict",
					species_name: 'Klebsiella pneumoniae'
				}
			},
		}
		'''
		# prevalence_sequence_id	model_id	species_name 		accession 	data_type		rgi_criteria	percent_identity
		# 10687						2882		Campylobacter coli	FBHN01		ncbi_contig		Strict			67.71
		# 0							1			2					3			4				5				6
		with open(self.model_species_data_type, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
			for row in reader:
				# add new model
				if row[1] not in variants.keys():
					variants.update({
						row[1]: { # model_id
							row[0]:  { # prevalence_sequence_id
								row[3]: #accession
									{
										"data_type":row[4],
										"rgi_criteria":row[5],
										"percent_identity":row[6],
										"species_name": row[2]
									}
							}
						}
					})
				# update existing model
				else:
					#  check if prev_id is present
					if row[0] not in variants[row[1]].keys():
						# new prevalence_sequence_id
						variants[row[1]].update({
							row[0]: { # prevalence_sequence_id
								row[3]: #accession
									{
										"data_type":row[4],
										"rgi_criteria":row[5],
										"percent_identity":row[6],
										"species_name": row[2]
									}
								}
						})
					else:
						# new accession
						variants[row[1]][row[0]].update({
								row[3]: #accession
									{
										"data_type":row[4],
										"rgi_criteria":row[5],
										"percent_identity":row[6],
										"species_name": row[2]
									}
						})

		return variants

	def get_baits_details(self):
		"""
		Parse index file to a dictionary for all baits
		"""
		baits = {}
		# 0			1		2		3		4		5		6			7
		# ProbeID, GeneID, TaxaID, ARO, ProbeSeq, Upstream, Downstream,RevComp
		with open(os.path.join(self.data, "baits-probes-with-sequence-info.txt"), 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter=',', quotechar='|')
			for row in reader:
				if row[0] != "ProbeID":
					baits.update({
						"{}|{}".format(row[0],row[3]): {
								"ProbeID": row[0],
								"GeneID":row[1],
								"TaxaID":row[2],
								"ARO": row[3],
								"ProbeSeq":row[4],
								"Upstream": row[5],
								"Downstream": row[6],
								"RevComp": row[7]
						}
					})
		return baits

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

	def filter_count_reads(self, is_mapped="false", length=""):
		"""
		Filter reads using mapQuality, length and isMapped
		TODO:: Review to remove / include
		"""
		read_one=os.path.join(self.working_directory, "{}.R1.fastq".format(self.output_file))
		read_two=os.path.join(self.working_directory, "{}.R2.fastq".format(self.output_file))
		options = ""

		if length:
			options = options + " -length >={}".format(length)

		if self.mapq:
			options = options + " -mapQuality >={}".format(self.mapq)

		filter_cmd = "bamtools filter -in {in_bam} -out {out_bam} -isMapped {is_mapped} {options}".format(
			in_bam=self.sorted_bam_sorted_file_length_100,
			out_bam=self.mapped,
			is_mapped=is_mapped,
			options=options
		)

		os.system(filter_cmd)

		# bedtools bamtobed -i trimmedreadstocardvslamq35l40.bam > trimmedreadstocardvslamq35l40.bed
		os.system("bedtools bamtobed -i {} > out.bed".format(self.mapped))

		extract_cmd = "samtools fastq -1 {read_one} -2 {read_two} {in_bam}".format(
			read_one=read_one,
			read_two=read_two,
			in_bam=self.mapped
		)
		os.system(extract_cmd)

		read_one_count_cmd = "awk '{c}' {read_one}".format(
			c="{s++}END{print s/4}",
			read_one=read_one
		)

		read_two_count_cmd = "awk '{c}' {read_two}".format(
			c="{s++}END{print s/4}",
			read_two=read_two
		)

		return os.popen(read_one_count_cmd).readlines()[0].strip("\n") ,  os.popen(read_two_count_cmd).readlines()[0].strip("\n")

	def get_stats(self):
		"""
		Get stats using 'bamtools stats' and 'samtools flagstat' function
		"""
		'''
		stats = {
			"mapped": {
				"read_one": 0,
				"read_two": 0
			},
			"unmapped": {
				"read_one": 0,
				"read_two": 0
			}
		}

		# unmapped
		unmapped = self.filter_count_reads()
		stats["unmapped"]["read_one"] = unmapped[0]
		stats["unmapped"]["read_two"] = unmapped[1]
		# mapped
		# max mapQuality for bowtie2 is 42
		# max mapQuality for bwa is 37
		#  see http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
		mapped = self.filter_count_reads(is_mapped="true")
		stats["mapped"]["read_one"] = mapped[0]
		stats["mapped"]["read_two"] = mapped[1]

		logger.info("count reads {}".format(json.dumps(stats,indent=2)))
		'''
		# overall stats for mapping
		cmd_overall = "bamtools stats -in {} > {}".format(self.sorted_bam_sorted_file_length_100 , self.mapping_overall_stats)
		# logger.info("overall mapping stats using {}".format(cmd_overall))
		# logger.info(cmd_overall)
		# os.system(cmd_overall)
		input_bam=self.sorted_bam_sorted_file_length_100
		output_stats=self.mapping_overall_stats
		subprocess.run(f"bamtools stats -in {input_bam} > {output_stats}", shell=True, check=True)
		# stats showing duplicates
		# write file reference_stats
		with open(self.mapping_artifacts_stats, "w") as out:
			out.write("**********************************************\n")
			out.write("Stats for Artifacts: \n")
			out.write("**********************************************\n")
			out.write("\n")
		cmd_artifacts = "samtools flagstat {} >> {}".format(self.sorted_bam_sorted_file_length_100 , self.mapping_artifacts_stats)
		# logger.info("mapping artifacts stats i.e duplicates using {}".format(cmd_artifacts))
		# logger.info(cmd_artifacts)
		# os.system(cmd_artifacts)
		input_bam=self.sorted_bam_sorted_file_length_100
		output_stats=self.mapping_artifacts_stats
		subprocess.run(f"samtools flagstat {input_bam} >> {output_stats}", shell=True, check=True)

	def find_between(s, start, end):
		return (s.split(start))[1].split(end)[0]

	def probes_stats(self, baits_card):
		stats = {}
		baits = {}
		with open(self.baits_mapping_data_tab, "r") as f2:
			reader=csv.reader(f2,delimiter='\t')
			for row in reader:
				if "ARO" in row[2]:
					bait = row[2]
					read = "{}|{}".format(row[0], row[1])

					if bait not in baits.keys():
						baits[bait] = [read]
					else:
						if read not in baits[bait]:
							baits[bait].append(read)

		aro_to_reads = {}

		for m in baits:
			aro = (m.split("|ARO:"))[1].split( "|")[0]
			if "|ARO:{}|".format(aro) in m:
				for r in baits[m]:
					if aro in aro_to_reads.keys():
						if r not in aro_to_reads[aro]:
							aro_to_reads[aro].append(r)
					else:
						aro_to_reads[aro] = [r]

		with open(self.aro_term_reads, "w") as tab_out3:
			writer = csv.writer(tab_out3, delimiter='\t', dialect='excel')
			writer.writerow([
				"ARO",
				"Number of Mapped Reads to baits"
				])
			for aro in aro_to_reads:
				writer.writerow([aro, len(aro_to_reads[aro])])

		with open(self.baits_mapping_data_json, "w") as outfile:
			json.dump(baits, outfile)

		reads_to_baits = {}
		for i in baits.keys():
			for j in baits[i]:
				t = i
				if j not in reads_to_baits.keys():
					reads_to_baits[j] = [t]
				else:
					if t not in reads_to_baits[j]:
					 reads_to_baits[j].append(t)

		with open(self.reads_mapping_data_json, "w") as outfile2:
			json.dump(reads_to_baits, outfile2)

		with open(self.baits_reads_count, "w") as tab_out2:
			writer = csv.writer(tab_out2, delimiter='\t', dialect='excel')
			writer.writerow([
				"Bait",
				"Number of Mapped Reads"
				])
			for item in baits:
				writer.writerow([item, len(baits[item])])

		probe_reads_count = {}

		with open(self.baits_reads_count, 'r') as csv_file:
			for row in csv.reader(csv_file, delimiter='\t'):
				probe_reads_count[row[0]] = row[1]

		data_out = {}
		for i in baits_card:
			data_out[i] = {}
			for k in baits_card[i]:
				probe_dict = k.split("|")
				probe = "|".join(probe_dict[:-1])
				if probe in probe_reads_count.keys():
					data_out[i].update({k : int(probe_reads_count[probe])})
				else:
					data_out[i].update({k : 0})

		with open(self.card_baits_reads_count_json, "w") as af:
			af.write(json.dumps(data_out,sort_keys=True))

		with open(self.reads_baits_count, "w") as tab_out2:
			writer = csv.writer(tab_out2, delimiter='\t', dialect='excel')
			writer.writerow([
				"Read",
				"Baits"
				])
			for item in reads_to_baits:
				writer.writerow([item,
				"; ".join(reads_to_baits[item])
				])

		with open(self.baits_mapping_data_tab, "r") as f2:
			reader=csv.reader(f2,delimiter='\t')
			for row in reader:
				if "ARO" in row[2]:
					term = row[2].split("|")[4]
					name = row[2].split("|")[5]
					probe = row[2]
					probe_ok = False
					for k in baits_card.keys():
						if term in k:
							matching = [s for s in baits_card[k] if probe in s]
							if matching:
								probe_ok = True

					if probe_ok == True:
						if term not in stats.keys():
							# read reference fasta to count baits used
							cmd = "cat {} | grep -c \"|{}|\"".format(self.reference_genome_baits,term)
							cmd2 = "cat {} | grep -c \"|{}|\"".format(self.reads_baits_count,term)
							stats[term] = {
								"aro_name": name,
								"total_baits": int(self.count_probes(cmd)),
								"read_count": int(self.count_probes(cmd2)),
								"mapped_baits": {
									probe: len(baits[probe])
								}
							}
						else:
							if probe not in stats[term]["mapped_baits"].keys():
								stats[term]["mapped_baits"].update({probe: len(baits[probe]) })

		# write tab for probes and mapped probes
		with open(self.mapping_baits_stats, "w") as tab_out:
			writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
			writer.writerow([
				"ARO Term",
				"ARO Accession",
				"Number of Baits",
				"Number of Mapped Baits with Reads",
				"Number of Reads Mapped to Baits",
				"Average Number of reads per Bait",
				"Number of reads per Bait Coefficient of Variation (%)"
				])

			'''
			Coefficient of variation = (Standard Deviation / Mean ) * 100
			- ratio of the standard devitation to the mean

			|    -     | Number of Probes | Mapped Probes |
			| -------- | --------         | ------------  |
			| Mean     |   57             |  22           |
			| std_dev  |   38             |  31           |
			| CV       | 66.67%           |  140.91%      |

			# note - probes with >0  mapping were considered

			'''

			for i in stats:
				accession, baits_count, baits_with_reads_count, reads_count, sample = self.get_counts(i, data_out)
				standard_devitation = 0
				coefficient_of_variation = 0

				try:
					standard_devitation = statistics.stdev(sample)
				except Exception as e:
					print(stats[i]["aro_name"], sample, e)

				mean = statistics.mean(sample)

				if mean > 0:
					coefficient_of_variation = (standard_devitation / mean) * 100

				writer.writerow([
					stats[i]["aro_name"],
					i,
					baits_count,
					baits_with_reads_count,
					reads_count,
					format(mean,'.2f'),
					format(coefficient_of_variation,'.2f')
				])


	def get_counts(self, accession, data_out):
		baits_count = 0
		baits_with_reads_count = 0
		reads_count = 0
		sample = []
		for i in data_out:
			if accession in i:
				baits_count = len(data_out[i])
				for c in data_out[i]:
					reads_count = reads_count + int(data_out[i][c])
					sample.append(int(data_out[i][c]))
					if int(data_out[i][c]) > 0:
						baits_with_reads_count = baits_with_reads_count + 1

				return accession, baits_count, baits_with_reads_count, reads_count, sample
		return accession, baits_count, baits_with_reads_count, reads_count, sample


	def count_probes(self, cmd):
		return os.popen(cmd).readlines()[0].strip("\n")

	def baits_reads_counts(self, accession):
		"""
		Returns
		number_of_mapped_baits,
		number_of_mapped_baits_with_reads,
		average_bait_coverage,
		bait_coverage_coefficient_of_variation
		"""
		if self.include_baits == True:
			with open(self.mapping_baits_stats, "r") as f2:
				reader=csv.reader(f2,delimiter='\t')
				for row in reader:
					if "ARO Term" not in row[0]:
						if accession in row[1]:
							return row[2], row[3], row[4], row[6]
				return 0, 0, 0, 0

		else:
			return 0, 0, 0, 0

	def get_model_id(self, models_by_accession, alignment_hit):
		model_id = ""
		if alignment_hit[0:22] == "Prevalence_Sequence_ID" or alignment_hit[0:4] == "ARO:":
			model_id = alignment_hit.split("|")[1].split(":")[1]
		else:
			accession = alignment_hit.split("|")[4].split(":")[1]
			try:
				model_id = models_by_accession[accession]["model_id"]
			except Exception as e:
				logger.warning("missing aro accession: {} for alignment {} -> {}".format(accession,alignment_hit,e))
		return model_id

	def get_mutation_details(self):
		# *.sam.temp.vcf.gz
		results = {}
		mutation_file_compressed = os.path.join(self.working_directory, "{}.temp.sam.temp.vcf.gz".format(self.output_file))
		# logger.debug("mutation_file_compressed: {}".format(mutation_file_compressed))

		with gzip.open(mutation_file_compressed, "rt") as m1:
			reader=csv.reader(m1,delimiter='\t')
			for row in reader:
				if "#" not in row[0]:
					if row[0] not in results.keys():
						results[row[0]] = []

					results[row[0]].append(
						"{reference_allele}{position}{alternative_allele}".format(
							reference_allele=row[3].strip(),
							position=row[1].strip(),
							alternative_allele=row[4].strip()
						)
					)

		return results

	def get_read_coverage(self):
		read_coverage_file = os.path.join(self.working_directory, "{}.temp.sam.temp.res".format(self.output_file))
		# logger.debug("read_coverage_file: {}".format(read_coverage_file))
		results = {}
		with open(read_coverage_file, "r") as c1:
				reader=csv.reader(c1,delimiter='\t')
				for row in reader:
					if row[0] != "#Template":
						results[row[0]] = { "depth": row[8].strip()}

		return results

	def get_consensus_sequence(self):
		consensus_sequence_file = os.path.join(self.working_directory, "{}.temp.sam.temp.fsa".format(self.output_file))
		# logger.debug("consensus_sequence_file: {}".format(consensus_sequence_file))
		result = []
		with open(consensus_sequence_file, "r") as c1:
			result = SeqIO.to_dict(SeqIO.parse(c1, "fasta"))
		return result


	def summary(self, alignment_hit, models, variants, baits, reads, models_by_accession,mutation,read_coverage,consensus_sequence):
		start = time.time()
		# logger.debug(alignment_hit)
		coverage = self.get_coverage_details(alignment_hit)
		model_id = self.get_model_id(models_by_accession, alignment_hit)

		try:
			alignments = self.get_alignments(alignment_hit)
			mapq_l = []
			mate_pair = []
			mapq_average = 0
			for a in alignments:
				mapq_l.append(int(a["mapq"]))
				if a["mrnm"] != "=" and a["mrnm"] not in mate_pair:
					if "ARO:{}".format(models[model_id]["ARO_accession"]) not in a["mrnm"]:
						mate_pair.append(a["mrnm"])

			if len(mapq_l) > 0:
				mapq_average = sum(mapq_l)/len(mapq_l)

			observed_in_genomes = "no data"
			observed_in_plasmids = "no data"
			prevalence_sequence_id = ""
			observed_data_types = []
			# range_of_reference_allele_source = []
			percent_identity = 0.0
			# Genus and species level only (only get first two words)
			observed_in_pathogens = []
			database = "CARD"
			reference_allele_source = "CARD curation"
			# logger.debug("alignment_hit: {}".format(alignment_hit))
			# if variants and "Resistomes & Variants" in database and "ARO:" not in alignment_hit:
			if "Prevalence_Sequence_ID" in alignment_hit:
				database = "Resistomes & Variants"
				# logger.debug("model_id: {}, alignment_hit: {}".format(model_id, alignment_hit))
				if model_id in variants.keys():
					_accession = ""
					for s in variants[model_id]:
						prevalence_sequence_id = alignment_hit.split("|")[0].split(":")[-1]
						observed_in_genomes = "NO"
						observed_in_plasmids = "NO"

						for accession in variants[model_id][prevalence_sequence_id]:
							_accession = accession

							if variants[model_id][prevalence_sequence_id][accession]["data_type"] not in observed_data_types:
								observed_data_types.append(variants[model_id][prevalence_sequence_id][accession]["data_type"])

							if variants[model_id][prevalence_sequence_id][accession]["species_name"] not in observed_in_pathogens:
								observed_in_pathogens.append(variants[model_id][prevalence_sequence_id][accession]["species_name"].replace('"', ""))

					if "Resistomes & Variants" in database:
						if "ncbi_chromosome" in observed_data_types:
							observed_in_genomes = "YES"
						if "ncbi_plasmid" in observed_data_types:
							observed_in_plasmids = "YES"

						# get prevalence_sequence_id Prevalence_Sequence_ID:10687|ID:2882|Name:tet(W/N/W)|ARO:3004442
						if "Prevalence_Sequence_ID" in alignment_hit:
							prevalence_sequence_id = alignment_hit.split("|")[0].split(":")[-1]
							try:
								reference_allele_source = "In silico {rgi_criteria} {percent_identity}% identity".format(
									rgi_criteria=variants[model_id][prevalence_sequence_id][_accession]["rgi_criteria"],
									percent_identity=variants[model_id][prevalence_sequence_id][_accession]["percent_identity"],
								)
								percent_identity = float(variants[model_id][prevalence_sequence_id][_accession]["percent_identity"])
							except Exception as e:
								reference_allele_source = ""
								# logger.debug(alignment_hit)
								# logger.debug(json.dumps(alignments, indent=2))
								# logger.debug(json.dumps(variants[model_id], indent=2))
								logger.warning("missing key with Prev_id: {}, Exception: {}, Database: {} for model_id: {}".format(prevalence_sequence_id, e, database, model_id))


				else:
					# provide info from model
					observed_in_pathogens = models[model_id]["taxon"]
					percent_identity = 100.0
			else:
				# logger.debug("model_id: {}, alignment_hit: {}".format(model_id, alignment_hit))
				observed_in_pathogens = models[model_id]["taxon"]
				# logger.debug(coverage)
				# assumption card canonical
				percent_identity = 100.0

			# check all clases categories
			resistomes = models[model_id]["categories"]
			if "AMR Gene Family" not in resistomes.keys():
				resistomes["AMR Gene Family"] = []
			if "Drug Class" not in resistomes.keys():
				resistomes["Drug Class"] = []
			if "Resistance Mechanism" not in resistomes.keys():
				resistomes["Resistance Mechanism"] = []

			stop = time.time()
			elapsed = stop - start
			logger.info("time lapsed: {} - {}".format(format(elapsed,'.3f'), alignment_hit))
			# self.async_print(alignment_hit, start, stop, elapsed)
			number_of_mapped_baits, number_of_mapped_baits_with_reads, average_bait_coverage, bait_coverage_coefficient_of_variation = self.baits_reads_counts(models[model_id]["ARO_accession"])

			consensus_sequence_dna = ""
			consensus_sequence_protein = ""
			read_coverage_depth = ""
			snps = ""

			if self.aligner == "kma":
				try:
					# if alignment_hit in consensus_sequence.keys():
					consensus_sequence_dna = "{}".format(consensus_sequence[alignment_hit].seq)
					input_seq = "{}".format(consensus_sequence[alignment_hit].seq)

					trailing_bases = len(input_seq) % 3

					if trailing_bases:
						input_seq = ''.join([input_seq, 'NN']) if trailing_bases == 1 else ''.join([input_seq, 'N'])

					try:
						consensus_sequence_protein = Seq.translate(input_seq,to_stop=False,table=11,cds=True)
					except Exception as e:
						consensus_sequence_protein = Seq.translate(input_seq,to_stop=False,table=11,cds=False)

					if trailing_bases:
						consensus_sequence_protein = consensus_sequence_protein[:-1]

					if alignment_hit in read_coverage.keys():
						read_coverage_depth = read_coverage[alignment_hit]["depth"]

					if alignment_hit in mutation.keys():
						snps = "; ".join(mutation[alignment_hit])

				except Exception as e:
					logger.warning("model with id : {}, has few mapped reads to make consensus sequence skipping: {}".format(model_id,e))
					return {}

			# report consensus if the depth is at least 5
			# logger.info("read_coverage_depth: {}".format(read_coverage_depth))
			if float(read_coverage_depth.strip() or 0) < 5.0:
				snps=""
				consensus_sequence_dna=""
				consensus_sequence_protein=""

			return {
				"id": alignment_hit,
				"cvterm_name": models[model_id]["model_name"],
				"aro_accession": models[model_id]["ARO_accession"],
				"model_type": models[model_id]["model_type"],
				"database": database,
				"reference_allele_source": reference_allele_source,
				"observed_in_genomes": observed_in_genomes,
				"observed_in_plasmids": observed_in_plasmids,
				"observed_in_pathogens": observed_in_pathogens,
				"range_of_reference_allele_source": percent_identity,
				"reads": reads[alignment_hit],
				"alignments": alignments,
				"mapq_average": format(mapq_average,'.2f'),
				"number_of_mapped_baits": number_of_mapped_baits,
				"number_of_mapped_baits_with_reads": number_of_mapped_baits_with_reads,
				"average_bait_coverage": average_bait_coverage,
				"bait_coverage_coefficient_of_variation": bait_coverage_coefficient_of_variation,
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
				,"predicted_pathogen": "N/A",
				"depth": read_coverage_depth, # get from *.sam.temp.res file when using kma
				"snps": snps,# get from vcf file
				"consensus_sequence_dna": consensus_sequence_dna, # get from *.sam.temp.fsa
				"consensus_sequence_protein": "{}".format(consensus_sequence_protein)
				}
		except Exception as e:
			logger.warning("missing model with id : {}, Exception: {}".format(model_id,e))

	def async_print(self, msg, start, stop, elapsed):
		logger.debug("{} ::: parent process: {} -> process id: {} ====|{}|{}|{}".format(
			msg, os.getppid(), os.getpid(),
			start,
			stop,
			elapsed
			)
		)

	def jobs(self, job):
		return self.summary(job[0], job[1], job[2], job[3], job[4], job[5], job[6], job[7], job[8])

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
		baits = {}
		models = {}
		mutation = {}
		read_coverage = {}
		consensus_sequence = {}
		logger.info("get_reads_count ...")
		reads = self.get_reads_count()
		logger.info("get_model_details ...")
		models, models_by_accession = self.get_model_details()

		if self.aligner == "kma":
			logger.info("get_mutation_details ...")
			mutation = self.get_mutation_details()
			logger.info("get_read_coverage ...")
			read_coverage = self.get_read_coverage()
			logger.info("get_consensus_sequence ...")
			consensus_sequence = self. get_consensus_sequence()

		if self.include_wildcard:
			logger.info("get_variant_details ...")
			variants = self.get_variant_details()
			# debug
			# with open("variants.json", "w") as outfile2:
			# 	json.dump(variants, outfile2)

		if self.include_baits:
			logger.info("get_baits_details ...")
			baits = self.get_baits_details()

		mapq_average = 0
		# t0 = time.time()

		jobs = []
		for alignment_hit in reads.keys():
			jobs.append((alignment_hit, models, variants, baits, reads, models_by_accession,mutation,read_coverage,consensus_sequence,))

		with Pool(processes=self.threads) as p:
			results = p.map_async(self.jobs, jobs)
			summary = results.get()

		# logger.info("Time: {}".format( format(time.time() - t0, '.3f')))
		# write json
		with open(self.allele_mapping_data_json, "w") as af:
			af.write(json.dumps(summary,sort_keys=True))

		# wrtie tab-delimited allele_mapping_data
		with open(self.allele_mapping_data_tab, "w") as tab_out:
			writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
			writer.writerow([
							"Reference Sequence",
							"ARO Term",
							"ARO Accession",
							"Reference Model Type",
							"Reference DB",
							"Reference Allele Source",
							"Resistomes & Variants: Observed in Genome(s)",
							"Resistomes & Variants: Observed in Plasmid(s)",
							"Resistomes & Variants: Observed Pathogen(s)",
							"Completely Mapped Reads",
							"Mapped Reads with Flanking Sequence",
							"All Mapped Reads",
							"Percent Coverage",
							"Length Coverage (bp)",
							"Average MAPQ (Completely Mapped Reads)",
							# "Number of Mapped Baits",
							# "Number of Mapped Baits with Reads",
							# "Average Number of reads per Bait",
							# "Number of reads per Bait Coefficient of Variation (%)",
							"Mate Pair Linkage",
							"Reference Length",
							# "Mutation",
							"AMR Gene Family",
							"Drug Class",
							"Resistance Mechanism",
							"Depth", # get from *.sam.temp.res file when using kma
							"SNPs", # get from vcf file
							"Consensus Sequence DNA", # get from *.sam.temp.fsa
							# ,"Predicted Pathogen"
							"Consensus Sequence Protein"
							])
			for r in summary:
				if r:
					writer.writerow([
						r["id"],
						r["cvterm_name"],
						r["aro_accession"],
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
						# r["number_of_mapped_baits"],
						# r["number_of_mapped_baits_with_reads"],
						# r["average_bait_coverage"],
						# r["bait_coverage_coefficient_of_variation"],
						"; ".join(r["mate_pair"]),
						r["reference"]["sequence_length"],
						# r["mutation"],
						"; ".join(r["resistomes"]["AMR Gene Family"]),
						"; ".join(r["resistomes"]["Drug Class"]),
						"; ".join(r["resistomes"]["Resistance Mechanism"])
						# ,r["predicted_pathogen"]
						,r["depth"]
						,r["snps"]
						,r["consensus_sequence_dna"]
						,r["consensus_sequence_protein"]
					])

		# wrtie tab-delimited gene_mapping_data
		mapping_summary = {}
		alleles_mapped = []
		index = "aro_accession"
		range_of_reference_allele_source = []
		for r in summary:
			if r:
				alleles_mapped.append(r[index])
				if r[index] not in mapping_summary.keys():
					mapping_summary[r[index]] = {
						"id": [],
						"cvterm_name": [],
						"aro_accession": [],
						"model_type": [],
						"database": [],
						"alleles_mapped": [],
						"range_of_reference_allele_source": [],
						"observed_in_genomes": [],
						"observed_in_plasmids": [],
						"observed_in_pathogens": [],
						"mapped": [],
						"unmapped": [],
						"all": [],
						"percent_coverage": [],
						"length_coverage": [],
						"mapq_average": [],
						"number_of_mapped_baits": [],
						"number_of_mapped_baits_with_reads": [],
						"average_bait_coverage": [],
						"bait_coverage_coefficient_of_variation": [],
						"mate_pair": [],
						"reference_sequence_length": [],
						"AMR Gene Family": [],
						"Drug Class": [],
						"Resistance Mechanism": []
					}

					mapping_summary[r[index]]["id"].append(r["id"])
					mapping_summary[r[index]]["cvterm_name"].append(r["cvterm_name"])
					mapping_summary[r[index]]["aro_accession"].append(r["aro_accession"])
					mapping_summary[r[index]]["model_type"].append(r["model_type"])
					mapping_summary[r[index]]["database"].append(r["database"])
					mapping_summary[r[index]]["observed_in_genomes"].append(r["observed_in_genomes"])
					mapping_summary[r[index]]["observed_in_plasmids"].append(r["observed_in_plasmids"])

					for p in r["observed_in_pathogens"]:
						mapping_summary[r[index]]["observed_in_pathogens"].append(p)

					mapping_summary[r[index]]["mapped"].append(r["reads"]["mapped"])
					mapping_summary[r[index]]["range_of_reference_allele_source"].append(r["range_of_reference_allele_source"])
					mapping_summary[r[index]]["unmapped"].append(r["reads"]["unmapped"])
					mapping_summary[r[index]]["all"].append(r["reads"]["all"])
					mapping_summary[r[index]]["percent_coverage"].append(r["percent_coverage"]["covered"])
					mapping_summary[r[index]]["length_coverage"].append(r["length_coverage"]["covered"])
					mapping_summary[r[index]]["mapq_average"].append(r["mapq_average"])
					mapping_summary[r[index]]["reference_sequence_length"].append(r["reference"]["sequence_length"])
					mapping_summary[r[index]]["number_of_mapped_baits"].append(r["number_of_mapped_baits"])
					mapping_summary[r[index]]["number_of_mapped_baits_with_reads"].append(r["number_of_mapped_baits_with_reads"])
					mapping_summary[r[index]]["average_bait_coverage"].append(r["average_bait_coverage"])
					mapping_summary[r[index]]["bait_coverage_coefficient_of_variation"].append(r["bait_coverage_coefficient_of_variation"])

					for m in r["mate_pair"]:
						if m not in ["*"]:
							arr = m.split("|")
							if len(arr) == 4:
								mapping_summary[r[index]]["mate_pair"].append("{}".format(m.split("|")[2].split(":")[1]))
							elif len(arr) == 7:
								mapping_summary[r[index]]["mate_pair"].append("{}".format(m.split("|")[5]))

					for a in r["resistomes"]["AMR Gene Family"]:
						mapping_summary[r[index]]["AMR Gene Family"].append(a)

					for d in r["resistomes"]["Drug Class"]:
						mapping_summary[r[index]]["Drug Class"].append(d)

					for c in r["resistomes"]["Resistance Mechanism"]:
						mapping_summary[r[index]]["Resistance Mechanism"].append(c)

				else:
					if r["model_type"] not in mapping_summary[r[index]]["model_type"]:
						mapping_summary[r[index]]["model_type"].append(r["model_type"])
					if r["database"] not in mapping_summary[r[index]]["database"]:
						mapping_summary[r[index]]["database"].append(r["database"])
					if r["observed_in_genomes"] not in mapping_summary[r[index]]["observed_in_genomes"]:
						mapping_summary[r[index]]["observed_in_genomes"].append(r["observed_in_genomes"])
					if r["observed_in_plasmids"] not in mapping_summary[r[index]]["observed_in_plasmids"]:
						mapping_summary[r[index]]["observed_in_plasmids"].append(r["observed_in_plasmids"])

					for p in r["observed_in_pathogens"]:
						if p not in mapping_summary[r[index]]["observed_in_pathogens"]:
							mapping_summary[r[index]]["observed_in_pathogens"].append(p)

					mapping_summary[r[index]]["mapped"].append(r["reads"]["mapped"])
					mapping_summary[r[index]]["range_of_reference_allele_source"].append(r["range_of_reference_allele_source"])
					mapping_summary[r[index]]["unmapped"].append(r["reads"]["unmapped"])
					mapping_summary[r[index]]["all"].append(r["reads"]["all"])
					mapping_summary[r[index]]["percent_coverage"].append(r["percent_coverage"]["covered"])
					mapping_summary[r[index]]["length_coverage"].append(r["length_coverage"]["covered"])
					mapping_summary[r[index]]["mapq_average"].append(r["mapq_average"])
					mapping_summary[r[index]]["reference_sequence_length"].append(r["reference"]["sequence_length"])
					mapping_summary[r[index]]["number_of_mapped_baits"].append(r["number_of_mapped_baits"])
					mapping_summary[r[index]]["number_of_mapped_baits_with_reads"].append(r["number_of_mapped_baits_with_reads"])
					mapping_summary[r[index]]["average_bait_coverage"].append(r["average_bait_coverage"])
					mapping_summary[r[index]]["bait_coverage_coefficient_of_variation"].append(r["bait_coverage_coefficient_of_variation"])

					for m in r["mate_pair"]:
						if m not in ["*"]:
							arr = m.split("|")
							if len(arr) == 4:
								mapping_summary[r[index]]["mate_pair"].append("{}".format(m.split("|")[2].split(":")[1]))
							elif len(arr) == 7:
								mapping_summary[r[index]]["mate_pair"].append("{}".format(m.split("|")[5]))

					for a in r["resistomes"]["AMR Gene Family"]:
						if a not in mapping_summary[r[index]]["AMR Gene Family"]:
							mapping_summary[r[index]]["AMR Gene Family"].append(a)

					for d in r["resistomes"]["Drug Class"]:
						if d not in mapping_summary[r[index]]["Drug Class"]:
							mapping_summary[r[index]]["Drug Class"].append(d)

					for c in r["resistomes"]["Resistance Mechanism"]:
						if c not in mapping_summary[r[index]]["Resistance Mechanism"]:
							mapping_summary[r[index]]["Resistance Mechanism"].append(c)

		with open(self.gene_mapping_data_tab, "w") as tab_out:
			writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
			writer.writerow([
							"ARO Term",
							"ARO Accession",
							"Reference Model Type",
							"Reference DB",
							"Alleles with Mapped Reads",
							"Reference Allele(s) Identity to CARD Reference Protein (%)",
							"Resistomes & Variants: Observed in Genome(s)",
							"Resistomes & Variants: Observed in Plasmid(s)",
							"Resistomes & Variants: Observed Pathogen(s)",
							"Completely Mapped Reads",
							"Mapped Reads with Flanking Sequence",
							"All Mapped Reads",
							"Average Percent Coverage",
							"Average Length Coverage (bp)",
							"Average MAPQ (Completely Mapped Reads)",
							"Number of Mapped Baits",
							"Number of Mapped Baits with Reads",
							"Average Number of reads per Bait",
							"Number of reads per Bait Coefficient of Variation (%)",
							"Number of reads mapping to baits and mapping to complete gene",
							"Number of reads mapping to baits and mapping to complete gene (%)",
							"Mate Pair Linkage (# reads)",
							"Reference Length",
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
						if k not in mapping_summary[i]["cvterm_name"]:
							mate_pairs.append("{} ({})".format(k,mp[k]))

				# identity range
				min_identity = float(min(mapping_summary[i]["range_of_reference_allele_source"]))
				max_identity = float(max(mapping_summary[i]["range_of_reference_allele_source"]))
				identity_range = ""
				# logger.debug("percent identity range for {} : {} => ({} - {})".format("; ".join(mapping_summary[i]["cvterm_name"]), mapping_summary[i]["range_of_reference_allele_source"], min_identity, max_identity))

				if min_identity == 0.0 and max_identity > 0.0:
					identity_range = "{}".format(max_identity)
				elif min_identity > 0.0 and max_identity > min_identity:
					identity_range = "{} - {}".format(min_identity, max_identity)
				elif min_identity == max_identity and min_identity > 0.0:
					identity_range = "{}".format(max_identity)

				# filter results
				if self.include_other_models is False:
					# show hits to homolog model only
					if "protein homolog model" in mapping_summary[i]["model_type"]:
						# logger.debug("show homolog only model types ...{}".format("; ".join(mapping_summary[i]["model_type"])))
						writer.writerow([
							"; ".join(mapping_summary[i]["cvterm_name"]),
							i,
							"; ".join(mapping_summary[i]["model_type"]),
							"; ".join(mapping_summary[i]["database"]),
							am[i],
							identity_range,
							observed_in_genomes,
							observed_in_plasmids,
							"; ".join(mapping_summary[i]["observed_in_pathogens"]),
							format(sum(map(float,mapping_summary[i]["mapped"])),'.2f'),
							format(sum(map(float,mapping_summary[i]["unmapped"])),'.2f'),
							format(sum(map(float,mapping_summary[i]["all"])),'.2f'),
							format(average_percent_coverage,'.2f'),
							format(average_length_coverage,'.2f'),
							format(average_mapq,'.2f'),
							mapping_summary[i]["number_of_mapped_baits"][-1],
							mapping_summary[i]["number_of_mapped_baits_with_reads"][-1],
							mapping_summary[i]["average_bait_coverage"][-1],
							mapping_summary[i]["bait_coverage_coefficient_of_variation"][-1],
							"N/A",
							"N/A",
							"; ".join(mate_pairs),
							"; ".join(mapping_summary[i]["reference_sequence_length"]),
							"; ".join(mapping_summary[i]["AMR Gene Family"]),
							"; ".join(mapping_summary[i]["Drug Class"]),
							"; ".join(mapping_summary[i]["Resistance Mechanism"])
						])

				else:
					# logger.debug("show all model types ... {}".format("; ".join(mapping_summary[i]["model_type"])))
					writer.writerow([
						"; ".join(mapping_summary[i]["cvterm_name"]),
						i,
						"; ".join(mapping_summary[i]["model_type"]),
						"; ".join(mapping_summary[i]["database"]),
						am[i],
						identity_range,
						observed_in_genomes,
						observed_in_plasmids,
						"; ".join(mapping_summary[i]["observed_in_pathogens"]),
						format(sum(map(float,mapping_summary[i]["mapped"])),'.2f'),
						format(sum(map(float,mapping_summary[i]["unmapped"])),'.2f'),
						format(sum(map(float,mapping_summary[i]["all"])),'.2f'),
						format(average_percent_coverage,'.2f'),
						format(average_length_coverage,'.2f'),
						format(average_mapq,'.2f'),
						mapping_summary[i]["number_of_mapped_baits"][-1],
						mapping_summary[i]["number_of_mapped_baits_with_reads"][-1],
						mapping_summary[i]["average_bait_coverage"][-1],
						mapping_summary[i]["bait_coverage_coefficient_of_variation"][-1],
						"N/A",
						"N/A",
						"; ".join(mate_pairs),
						"; ".join(mapping_summary[i]["reference_sequence_length"]),
						"; ".join(mapping_summary[i]["AMR Gene Family"]),
						"; ".join(mapping_summary[i]["Drug Class"]),
						"; ".join(mapping_summary[i]["Resistance Mechanism"])
					])

	def check_index(self, index_directory, reference_genome):
		"""
		Check if index exists for a given reference fasta file.
		"""
		logger.info("check database index")
		if self.aligner == "bowtie2":
			files = [os.path.basename(x) for x in glob.glob(os.path.join(os.path.dirname(index_directory),"*"))]
			logger.info(json.dumps(files, indent=2))
			if (("bowtie2.1.bt2" in files) and \
				("bowtie2.2.bt2" in files) and \
				("bowtie2.3.bt2" in files) and \
				("bowtie2.4.bt2" in files) and \
				("bowtie2.rev.1.bt2" in files) and \
				("bowtie2.rev.2.bt2" in files)) == False:
				# create index and save results in ./db from reference genome: (.fasta)
				logger.info("create index for reference: {} using aligner: {} ".format(reference_genome, self.aligner))
				self.create_index(index_directory=index_directory,reference_genome=reference_genome)
			else:
				logger.info("index already exists for reference: {} using aligner: {}".format(reference_genome,self.aligner))
		elif self.aligner == "kma":
			files = [os.path.basename(x) for x in glob.glob(os.path.join(os.path.dirname(index_directory),"*"))]
			logger.info(json.dumps(files, indent=2))
			if (("kma.comp.b" in files) and \
				("kma.length.b" in files) and \
				("kma.seq.b" in files) and \
				("kma.name" in files)) == False:
				# create index and save results in ./db from reference genome: (.fasta)
				logger.info("create index for reference: {} using aligner: {} ".format(reference_genome, self.aligner))
				self.create_index(index_directory=index_directory,reference_genome=reference_genome)
			else:
				logger.info("index already exists for reference: {} using aligner: {}".format(reference_genome,self.aligner))
		else:
			files = [os.path.basename(x) for x in glob.glob(os.path.join(os.path.dirname(index_directory),"*"))]
			logger.info(json.dumps(files, indent=2))
			if (("bwa.amb" in  files) and \
				("bwa.ann" in  files) and \
				("bwa.bwt" in  files) and \
				("bwa.pac" in  files) and \
				("bwa.sa" in files)) == False:
				# create index and save results in ./db from reference genome: (.fasta)
				logger.info("create index for reference: {} using aligner: {} ".format(reference_genome, self.aligner))
				self.create_index(index_directory=index_directory,reference_genome=reference_genome)
			else:
				logger.info("index already exists for reference: {} using aligner: {}".format(reference_genome, self.aligner))

	def run(self):
		"""
		Align reads to reference genomes and report
		"""
		t0 = time.time()
		logger.info("inputs")
		logger.info(json.dumps(self.__dict__, indent=2))

	    # check index / create index / align
		logger.info("align using {}".format(self.aligner))

		if self.aligner == "bowtie2":
			if self.read_two == None:
				self.align_bowtie2_unpaired(reference_genome=self.reference_genome, index_directory=self.index_directory_bowtie2, output_sam_file=self.output_sam_file)
			else:
				self.align_bowtie2(reference_genome=self.reference_genome, index_directory=self.index_directory_bowtie2, output_sam_file=self.output_sam_file)
		elif self.aligner == "kma":
			if self.read_two == None:
				self.align_kma_interleaved(reference_genome=self.reference_genome, index_directory=self.index_directory_kma, output_sam_file=self.output_sam_file)
			else:
				self.align_kma(reference_genome=self.reference_genome, index_directory=self.index_directory_kma, output_sam_file=self.output_sam_file)
		else:
			if self.read_two == None:
				self.align_bwa_single_end_mapping(reference_genome=self.reference_genome, index_directory=self.index_directory_bwa, output_sam_file=self.output_sam_file)
			else:
				self.align_bwa_paired_end_mapping(reference_genome=self.reference_genome, index_directory=self.index_directory_bwa,  output_sam_file=self.output_sam_file)

		# convert SAM file to BAM file
		logger.info("convert SAM file to BAM file")
		self.convert_sam_to_bam(input_sam_file=self.output_sam_file, output_bam_file=self.output_bam_file)

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

		# pull aligned
		logger.info("pull aligned")
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

		if self.include_baits == True:

			# map baits to complete genes
			logger.debug("map baits to complete genes")
			self.align_bowtie2_baits_to_genes(
				reference_genome=self.reference_genome,
				index_directory=self.index_directory_bowtie2,
				output_sam_file=self.baits_card_sam
			)

			logger.info("convert SAM file to BAM file")
			self.convert_sam_to_bam(input_sam_file=self.baits_card_sam, output_bam_file=self.baits_card_bam)

			os.system("samtools view -F4 --threads {threads} {input_bam} | cut -f 1,2,3 | sort -s -n -k 1,1 > {output_tab}".format(
				threads=self.threads,
				input_bam=self.baits_card_bam,
				output_tab=self.baits_card_data_tab
				))

			baits_card = {}
			with open(self.baits_card_data_tab, "r") as f2:
				reader=csv.reader(f2,delimiter='\t')
				for row in reader:
					gene = row[2]
					bait = "{}|{}".format(row[0], row[1])
					if gene not in baits_card.keys():
						baits_card[gene] = [bait]
					else:
						baits_card[gene].append(bait)

			# write json
			with open(self.baits_card_json, "w") as af:
				af.write(json.dumps(baits_card,sort_keys=True))

			with open(self.baits_card_tab, "w") as tab_out:
				writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
				writer.writerow(["Gene", "Number of baits mapped to gene"])
				for g in baits_card:
					writer.writerow([g,len(baits_card[g])])

			# map reads to baits
			logger.debug("map reads to baits...")
			self.align_bowtie2(reference_genome=self.reference_genome_baits, index_directory=self.index_directory_bowtie2_baits, output_sam_file=self.output_sam_file_baits)

			logger.info("convert SAM file to BAM file")
			self.convert_sam_to_bam(input_sam_file=self.output_sam_file_baits, output_bam_file=self.output_bam_file_baits)

			# get mapped
			logger.info("get number of reads mapped to baits")
			os.system("samtools view -F4 --threads {threads} {input_bam} | cut -f 1,2,3 | sort -s -n -k 1,1 > {output_tab}".format(
				threads=self.threads,
				input_bam=self.output_bam_file_baits,
				output_tab=self.baits_mapping_data_tab
				))

			# get stats
			logger.debug("get baits statistics")
			self.probes_stats(baits_card=baits_card)

		# get summary
		logger.info("get summary")
		self.get_summary()

		# get stats
		logger.info("get statistics")
		self.get_stats()

		# clean temporary files
		logger.info("clean temporary files")
		self.clean_files()
		logger.info("Total Time: {} seconds".format( format(time.time() - t0, '.3f')))
		logger.info("Done.")
