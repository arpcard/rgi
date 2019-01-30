from app.settings import os, SeqIO, logger
import tempfile, time, fileinput, math, multiprocessing, shutil

class ORF(object):
	"""Class to find open reading frames from nucleotide sequence."""
	def __init__(self,input_file, threads, clean=True, working_directory=None, low_quality=False, training_file=None, split_prodigal_jobs=False):
		"""Creates ORF object for finding open reading frames."""
		self.input_file = input_file
		self.clean = clean
		self.working_directory = working_directory
		self.low_quality = low_quality
		self.training_file = training_file
		self.threads = threads
		self.split_prodigal_jobs = split_prodigal_jobs

	def __repr__(self):
		"""Returns ORF class full object."""
		return "ORF({}".format(self.__dict__)

	def contig_to_orf(self):
		"""Converts contigs to open reading frames."""
		if self.training_file != None:
			self.orf_prodigal_train()
		else:
			self.orf_prodigal()

	def min_max_sequence_length(self):
		"""Returns minimum and maximun sequence length in multi-fasta inputs""" 
		sequences = []
		for record in SeqIO.parse(self.input_file, "fasta"):
			sequences.append(len(record.seq))
		return min(sequences), max(sequences), len(sequences)

	def orf_prodigal(self):
		"""Runs PRODIGAL to find open reading frames."""
		quality = "-n -p single"

		minimum_sequence_length, maximum_sequence_length, number_of_sequences = self.min_max_sequence_length()
		logger.info("minimum sequence length: {}, maximun sequence length {}, number of sequences: {}".format(minimum_sequence_length,maximum_sequence_length, number_of_sequences))

		if number_of_sequences > 1 and self.split_prodigal_jobs == True:
			# TODO validate if fasta file doesn't contain gaps
			self.orf_prodigal_multi()
		else:
			if self.low_quality == True or minimum_sequence_length < 20000:
				quality = "-p meta"

			filename = os.path.basename(self.input_file)

			stdout = "2> /dev/null"

			cmd = "prodigal -q -m -a {trans_file} -i {input_file} -o  {output_file} -d {nuc_file} -s {potential_genes} {quality} {stdout}" \
			.format(
					trans_file=os.path.join(self.working_directory, "{}.temp.contig.fsa".format(filename)),
					input_file=self.input_file,
					output_file=os.path.join(self.working_directory, "{}.temp.draft".format(filename)),
					quality=quality,
					stdout=stdout,
					nuc_file= os.path.join(self.working_directory,  "{}.temp.contigToORF.fsa".format(filename)),
					potential_genes= os.path.join(self.working_directory,  "{}.temp.potentialGenes".format(filename))
				)

			# logger.debug(cmd)
			os.system(cmd)

			# format the contig file headers to remove space
			#format_fasta_headers(working_directory+"/"+filename+".contig.fsa")

			if self.clean == True:
				os.remove(os.path.join(self.working_directory, "{}.temp.draft".format(filename)))

	def orf_prodigal_multi(self):
		seq = self.split_fasta()
		self.execute_threads(seq)

	def worker(self, input_fasta):
		o_f_path, o_f_name = os.path.split(os.path.abspath(input_fasta))
		cmd = "prodigal -p meta -q -m -i {input_fasta} -d {wd}/{tmp_name}.temp.contigToORF.fsa \
		-a {wd}/{tmp_name}.temp.contig.fsa \
		-o {wd}/{tmp_name}.temp.draft \
		-s {wd}/{tmp_name}.temp.potentialGenes 2> /dev/null".format(input_fasta=input_fasta,wd=self.working_directory,tmp_name=o_f_name)
		# logger.debug(cmd)
		os.system(cmd)

	def prodigal_run(self, fasta, *o):
		files = []
		logger.info("prodigal_run on {} sequences...".format(len(fasta)))
		o_f_path, o_f_name = os.path.split(os.path.abspath(self.input_file))
		for entry in fasta:
			# create directory tmp if it doesn't exist
			tmp = os.path.join(self.working_directory, "{}.temp.directory".format(o_f_name))
			if not os.path.exists(tmp):
				os.makedirs(tmp, exist_ok=True)
			with tempfile.NamedTemporaryFile(mode='w+', dir=tmp, delete=False) as fp:
				fp.write(">{}\n{}\n".format(entry.id, entry.seq))
				if fp.closed == False:
					fp.close()
					files.append(fp.name)
					self.worker(fp.name)
		o[0].put(files)

	def write_output_file(self, output, filenames):
		with open(os.path.join(self.working_directory, output), 'a+') as fout, fileinput.input(filenames) as fin:
			for line in fin:
				fout.write(line)
		if os.path.exists(filenames[0]):
			if self.clean == True:
				logger.debug("Removed temp file: {}".format(filenames[0]))
				os.remove(filenames[0])

	def chunk_list(self, iterator, n):
		"""
		Reference: https://biopython.org/wiki/Split_large_file
		"""
		entry = True
		while entry:
			batch = []
			while len(batch) < n:
				try:
					entry = next(iterator)
				except StopIteration:
					entry = None
				if entry is None: # end of file
					break
				batch.append(entry)
			if batch:
				yield batch

	def split_fasta(self):
		iterator = SeqIO.parse(self.input_file, "fasta")
		temp_iterator = SeqIO.parse(self.input_file, "fasta")
		# counts number of sequences
		ns = sum(1 for i in temp_iterator)
		# maximizes list size for threads available
		list_size = math.ceil(ns/self.threads)
		# returns a list of lists
		split_sequences = list(self.chunk_list(iterator, list_size))
		return split_sequences

	def execute_threads(self, split_sequences):
		output = multiprocessing.Queue()
		processes = []
		for ind in range(len(split_sequences)):
			process = multiprocessing.Process(target=self.prodigal_run,args=(split_sequences[ind], output,))
			process.start()
			processes.append(process)
		results = [output.get() for process in processes]
		for process in processes:
			process.join()

		f_path, f_name = os.path.split(os.path.abspath(self.input_file))

		output_dna_orf = "{wd}/{tmp_name}.temp.contigToORF.fsa".format(wd=self.working_directory, tmp_name=f_name)
		output_prot_orf = "{wd}/{tmp_name}.temp.contig.fsa".format(wd=self.working_directory,tmp_name=f_name)
		output_draft = "{wd}/{tmp_name}.temp.draft".format(wd=self.working_directory,tmp_name=f_name)
		output_potential_genes = "{wd}/{tmp_name}.temp.potentialGenes".format(wd=self.working_directory,tmp_name=f_name)

		# combine results
		for i in results:
			for j in i:
				o_f_path, o_f_name = os.path.split(os.path.abspath(j))
				self.write_output_file(output_dna_orf, ["{wd}/{tmp_name}.temp.contigToORF.fsa".format(wd=self.working_directory,tmp_name=o_f_name)])
				self.write_output_file(output_prot_orf, ["{wd}/{tmp_name}.temp.contig.fsa".format(wd=self.working_directory,tmp_name=o_f_name)])
				self.write_output_file(output_draft, ["{wd}/{tmp_name}.temp.draft".format(wd=self.working_directory,tmp_name=o_f_name)])
				self.write_output_file(output_potential_genes, ["{wd}/{tmp_name}.temp.potentialGenes".format(wd=self.working_directory,tmp_name=o_f_name)])
		# remove temps directories
		tmp = os.path.join(self.working_directory, "{}.temp.directory".format(f_name))
		if os.path.exists(tmp):
			if self.clean == True:
				logger.debug("Removed directory: {}".format(tmp))
				shutil.rmtree(tmp)
		return results

	def orf_prodigal_train(self):
		"""Runs PRODIGAL to find open reading frames using a training file from complete genomes references."""
		training_file = os.path.join(self.training_file)

		if os.path.exists(training_file):
			quality = " -t {} ".format(training_file)
			filename = os.path.basename(self.input_file)
			stdout = "2> /dev/null"
			cmd = "prodigal -q -m -a {trans_file} -i {input_file} -o  {output_file} -d {nuc_file} -s {potential_genes} {quality} {stdout}" \
			   .format(
					trans_file=os.path.join(self.working_directory, "{}.temp.contig.fsa".format(filename)),
					input_file=self.input_file,
					output_file=os.path.join(self.working_directory, "{}.temp.draft".format(filename)),
					quality=quality,
					stdout=stdout,
					nuc_file= os.path.join(self.working_directory,  "{}.temp.contigToORF.fsa".format(filename)),
					potential_genes= os.path.join(self.working_directory,  "{}.temp.potentialGenes".format(filename))
				)
			# logger.debug(cmd)
			os.system(cmd)

			if self.clean == True:
				os.remove(os.path.join(self.working_directory, "{}.temp.draft".format(filename)))
		else:
			logger.error("Missing training file: {} ".format(training_file))
			exit()

	def get_character_len(self,file_path):
		"""Returns character count in a file."""
		chars = words = lines = 0
		with open(file_path, 'r') as in_file:
		    for line in in_file:
		    	if line[0] == '>':
		    		pass
		    	else:
			        lines += 1
			        words += len(line.split())
			        chars += len(line)
		# logger.info("chars count: {}".format(chars))
		return chars




