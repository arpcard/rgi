from app.settings import *

class Blast(object):
	"""Class to create Blast object and align for protein and translated DNA searches."""
	def __init__(self,input_file, output_file=None, program = 'blastp', num_threads=4 ):
		"""Creates Blast object for running NCBI BLAST algorithm."""
		self.input_file = input_file
		if output_file == None:
			f_path, f_name = os.path.split(input_file)
			self.output_file = os.path.join(f_path,"{}.blastRes.xml".format(f_name))
		else:
			self.output_file = output_file
		self.program = program
		self.num_threads = num_threads
		self.outfmt = 5

	def __repr__(self):
		"""Returns Blast class full object."""
		return "Blast({}".format(self.__dict__)

	def run(self):
		"""Runs BLAST algorithm.""" 
		logger.info("run blast")
		os.system('{program} -query {input} -db {path}protein.db \
					-num_threads {num_threads} -outfmt {outfmt} -out {output_file}' \
					.format(
						program=self.program, 
						num_threads=self.num_threads, 
						outfmt=self.outfmt,
						input=self.input_file,
						path=path,
						output_file=self.output_file
					)
				)

	def run_custom(self, db ):
		"""Runs DIAMOND algorithm.""" 
		# logger.info("run diamond")
		os.system('{program} -query {input} -db {db} -num_threads {num_threads} \
					-outfmt {outfmt} -out {output_file} -perc_identity {perc_identity} \
					-strand {strand}' \
					.format(
						program=self.program, 
						input=self.input_file,
						db=db,
						num_threads=self.num_threads, 
						outfmt=self.outfmt,
						output_file=self.output_file,
						perc_identity=95.0,
						strand="both"
					)
				)
		logger.info("done running {} -> {}".format(self.program, db))


