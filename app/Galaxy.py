from app.settings import *

class Galaxy(object):
	def __init__(self, database):
		self.database = database

	def load_db_galaxy(self):
		card_dir = self.database

		logger.info("source_path: ", card_dir)

		if card_dir == None:
			logger.error("no new card path")
			exit("no new card path")

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

		files = os.listdir(card_dir)

		for f in files:
			if not f.startswith('.'):
				found_files.append(f)

		missing_files = list(set(needed_files) - set(found_files))

		if len(missing_files) > 0 :
		 	logger.error("Error: missing database files {}".format(missing_files))
		 	exit("Error: missing database files {}".format(missing_files))
		
		# Files found - move files into _data and _db directories
		for f in os.listdir(card_dir):
			if not f.startswith('.') and f in needed_files:
				src_path = os.path.join(card_dir,f)
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

				
