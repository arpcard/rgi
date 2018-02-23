from app.Base import RGIBase
from app.Database import Database
from app.Blast import Blast
from app.Diamond import Diamond
from app.ORF import ORF
from app.Filter import Filter

import filetype
from Bio import SeqIO
import glob
import time

from app.settings import *

class RGI(RGIBase):
	"""Class to predict resistome(s) from protein or nucleotide data based on CARD detection models."""

	def __init__(self,input_type='contig',input_sequence=None,threads=32,output_file=None,loose=False, \
				clean=True,data='na',aligner='blast',galaxy=None, local_database=False, low_quality=False):
		"""Creates RGI object for resistome(s) prediction."""

		o_f_path, o_f_name = os.path.split(os.path.abspath(output_file))

		self.input_type = input_type.lower()
		self.input_sequence = os.path.abspath(input_sequence)
		self.threads = threads
		self.num_sequences = 1
		self.output_file = os.path.abspath(output_file)
		self.loose = loose
		self.clean = clean
		self.data = data
		self.aligner = aligner.lower()
		self.database = galaxy
		self.low_quality = low_quality

		self.local_database = local_database
		self.db = path
		self.dp = data_path

		if self.local_database:
			self.db = LOCAL_DATABASE
			self.dp = LOCAL_DATABASE

		self.working_directory = o_f_path
		self.blast_results_xml_file = ''

		super(RGIBase, self).__init__()

	def __repr__(self):
		"""Returns RGI class full object."""
		return "RGI({}".format(self.__dict__)

	@classmethod
	def from_string(cls,cmd_string):
		"""Creates RGI object from string."""
		input_type,input_sequence,threads,num_sequences,output_file,aligner,database = cmd_string.split('@')
		return cls(input_type,input_sequence,threads,num_sequences,output_file,aligner,database)

	@classmethod
	def from_args(cls, *initial_data, **kwargs):
		"""Creates RGI object from args."""
		for dictionary in initial_data:
			for key in dictionary:
				if key in ['input_type','loose','clean','aligner']:
					setattr(cls, key, dictionary[key].lower())
				setattr(cls, key, dictionary[key])

		for key in kwargs:
			if key in ['input_type','loose','clean','aligner']:
				setattr(cls, key, kwargs[key].lower())
			setattr(cls, key, kwargs[key])

		return cls()

	def validate_inputs(self):
		"""Validate inputs."""
		# give warning if input_type is neither 'protein' nor 'dna' and then exit this program
		if self.input_type != 'protein' and self.input_type != 'contig' and self.input_type != 'read':
			logger.error("input_type must be one of protein, contig or read.")
			exit()
		if self.input_sequence == None:
			logger.error("missing input file")
			exit()

		logger.info("{} => {}".format(self.input_sequence, filetype.guess(self.input_sequence)))
		kind = filetype.guess(self.input_sequence)

		if kind is None:
			if self.is_fasta(self.input_sequence):
				logger.info("Fasta")
		else:
			logger.error(kind.extension)
			logger.error(kind.mime)
			logger.warning("Sorry, no supoprt for this format.")
			exit()

	@staticmethod
	def is_fasta(filename):
	    with open(filename, "r") as handle:
	        fasta = SeqIO.parse(handle, "fasta")
	        return any(fasta)

	
	def __set_xml_filepath(self,fp):
		"""Sets blast xml filepath."""
		self.blast_results_xml_file = fp

	def create_databases(self):
		"""Creates databases."""
		db_obj = Database(self.local_database)
		db_obj.build_databases()

	def run(self):
		"""Runs RGI."""
		t0 = time.time()
		self.validate_inputs()
		self.create_databases()
		self.run_blast()
		self.filter_process()
		#logger.info("Output......")
		#self.out()
		logger.info('Total running time {}s'.format(round(time.time() - t0, 3)))

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
			if os.path.basename(self.input_sequence) + ".temp" in f and os.path.isfile(f):
				self.remove_file(f)				

	def remove_file(self, f):
		"""Removes file."""
		if os.path.exists(f):
			logger.info("Removed file: {}".format(f))
			os.remove(f)
		else:
			logger.warning("Missing file: {}".format(f))

	def out(self):
		"""Writes tab-delimited, ggf3 output files."""
		tab_obj = Output(self.output_file)
		tab_obj.run()

	def run_blast(self):
		"""Runs blast."""
		if self.input_type == "protein":
			self.process_protein()
		elif  self.input_type == "contig":
			self.process_contig()
		elif self.input_type == "read":
			self.process_read()
			# logger.debug("TODO:: add read functions")
			exit("")
		else:
			logger.error("Invalid input_type: {} ".format(self.input_type))
			exit()

	def set_xml_filepath(self,fp):
		"""Sets blast xml filepath."""
		logger.info("set blast xml file: [{}]".format(fp))
		self.blast_results_xml_file = fp

	def process_protein(self):
		"""Process protein sequence(s)."""
		file_name = os.path.basename(self.input_sequence)
		xml_file = os.path.join(self.working_directory,"{}.temp.blastRes.xml".format(file_name))

		if self.aligner == "diamond":
			diamond_obj = Diamond(self.input_sequence, xml_file, local_database=self.local_database, num_threads=self.threads)
			diamond_obj.run()
		else:
			blast_obj = Blast(self.input_sequence, xml_file, local_database=self.local_database, num_threads=self.threads)
			blast_obj.run()

		self.set_xml_filepath(xml_file)

	def process_contig(self):
		"""Process nuclotide sequence(s)."""
		file_name = os.path.basename(self.input_sequence)
		orf_obj = ORF(input_file=self.input_sequence, clean=self.clean, working_directory=self.working_directory, low_quality=self.low_quality)
		orf_obj.contig_to_orf()
		contig_fsa_file = os.path.join(self.working_directory,"{}.temp.contig.fsa".format(file_name))
		blast_results_xml_file = os.path.join(self.working_directory,"{}.temp.contig.fsa.blastRes.xml".format(file_name))

		try:
			if os.stat(contig_fsa_file).st_size > 0:
				logger.info("work with file {}".format(contig_fsa_file))
				if self.aligner == "diamond":
					diamond_obj = Diamond(contig_fsa_file, local_database=self.local_database, num_threads=self.threads)
					diamond_obj.run()
				else:
					blast_obj = Blast(contig_fsa_file, local_database=self.local_database, num_threads=self.threads)
					blast_obj.run()
				self.set_xml_filepath(blast_results_xml_file)
			else:
				self.write_stub_output_file()
				logger.error("no open reading frames (orfs) found.")
		except Exception as e:
			self.write_stub_output_file()
			logger.exception("failed to write orf file")
		else:
			# logger.info("success procession orf file")
			pass

	def write_stub_output_file(self):
		# write empty output file if there are no open reading frames
		with open(os.path.join(self.output_file), 'w') as fout:
			fout.write(json.dumps({}))

	def process_read(self):
		"""Process fastq or reads sequence(s)."""
		logger.info("process read")

	# @profile
	def filter_process(self):
		logger.info("run filter")
		"""Filter each detection models and predict resistome(s)."""
		filter_obj = Filter(self.input_type,  self.loose, self.input_sequence, self.blast_results_xml_file, \
			os.path.join(self.dp,"card.json"),os.path.basename(self.input_sequence) ,self.output_file, self)
		filter_obj.run()

	def output(self): pass




			

