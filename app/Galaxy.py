from app.settings import *

class Galaxy(object):
	def __init__(self, database, debug=False):
		self.database = database
		self.debug = debug
		if self.debug:
			logger.setLevel(10)

	def load_db_galaxy(self):
		logger.info("source_path: {}".format(self.database))

		if self.database == None:
			logger.error("no new card path")
			exit()

		"""
		verify that we have the following files at the specified location:

		# CARD json data file
		- card.json

		# diamond blast
		- protein.db.dmnd

		# ncbi blast
		- protein.db.phr
		- protein.db.pin
		- protein.db.psq

		# Protein fasta file
		- proteindb.fsa

		"""
		needed_files = ['card.json','proteindb.fsa','protein.db.dmnd','protein.db.phr','protein.db.pin','protein.db.psq']
		found_files = []

		files = os.listdir(self.database)

		for f in files:
			if not f.startswith('.'):
				found_files.append(f)

		missing_files = list(set(needed_files) - set(found_files))

		if len(missing_files) > 0:
		 	logger.error("Error: missing database files {}".format(missing_files))
		 	exit()
		
		# Files found - move files into _data and _db directories
		for f in os.listdir(self.database):
			if not f.startswith('.') and f in needed_files:
				src_path = os.path.join(self.database,f)
				dst_path = ""
				if f in ['card.json']:
					dst_path = os.path.join(data_path,f)
				else:
					dst_path = os.path.join(path,f)
				logger.info("copy {} to {}".format(src_path,dst_path))
				try:
					shutil.copy2(src_path,dst_path)
				except Exception as e:
					logger.warning("failed to copy file: {}".format(e))

				
