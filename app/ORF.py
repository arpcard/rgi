from app.settings import os

class ORF(object):
	"""Class to find open reading frames from nucleotide sequence."""
	def __init__(self,input_file, clean=True, working_directory=None):
		"""Creates ORF object for finding open reading frames."""
		self.input_file = input_file
		self.clean = clean
		self.working_directory = working_directory

	def __repr__(self):
		"""Returns ORF class full object."""
		return "ORF({}".format(self.__dict__)

	def contig_to_orf(self):
		"""Converts contigs to open reading frames."""
		self.orf_prodigal()

	def orf_prodigal(self):
		"""Runs PRODIGAL to find open reading frames."""
		# logger.info("get open reading frames using prodigal")
		p = "single"
		count =  int(self.get_character_len(self.input_file))
		#check for number of characters for selection procedure -- ribosomal binding sites
		if count < 200000:
			p = "meta"

		filename = os.path.basename(self.input_file)
		os.system("prodigal -c -m -a {trans_file}  -i {input_file} -o  {output_file} -p {mode} -d {nuc_file} -q -s {potential_genes}" \
		   .format(
				trans_file=os.path.join(self.working_directory, "{}.temp.contig.fsa".format(filename)),
				input_file=self.input_file,
				output_file=os.path.join(self.working_directory, "{}.temp.draft".format(filename)),
				mode=p,
				nuc_file= os.path.join(self.working_directory,  "{}.temp.contigToORF.fsa".format(filename)),
				potential_genes= os.path.join(self.working_directory,  "{}.temp.potentialGenes".format(filename))
			)
		)
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




