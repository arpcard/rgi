from app.settings import os, SeqIO, logger

class ORF(object):
	"""Class to find open reading frames from nucleotide sequence."""
	def __init__(self,input_file, clean=True, working_directory=None, low_quality=False):
		"""Creates ORF object for finding open reading frames."""
		self.input_file = input_file
		self.clean = clean
		self.working_directory = working_directory
		self.low_quality = low_quality

	def __repr__(self):
		"""Returns ORF class full object."""
		return "ORF({}".format(self.__dict__)

	def contig_to_orf(self):
		"""Converts contigs to open reading frames."""
		self.orf_prodigal()

	def min_max_sequence_length(self):
		"""Returns minimum and maximun sequence length in multi-fasta inputs""" 
		sequences = []
		for record in SeqIO.parse(self.input_file, "fasta"):
			sequences.append(len(record.seq))
		return min(sequences), max(sequences)

	def orf_prodigal(self):
		"""Runs PRODIGAL to find open reading frames."""
		quality = "-n -p single"

		_min, _max = self.min_max_sequence_length()
		logger.info("minimum sequence length: {}, maximun sequence length {}".format(_min,_max))

		if self.low_quality == True or _min < 20000:
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

		logger.info(cmd)
		os.system(cmd)

		# format the contig file headers to remove space
		#format_fasta_headers(working_directory+"/"+filename+".contig.fsa")

		if self.clean == True:
			os.remove(os.path.join(self.working_directory, "{}.temp.draft".format(filename)))


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




