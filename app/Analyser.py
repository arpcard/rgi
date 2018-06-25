from app.settings import *

class Analyser(object):
	"""Class to find analyse fasta files."""
	def __init__(self,input_file):
		"""Creates Analyser object for analysing fasta files."""
		self.input_file = input_file
		logger.info(repr(self))

	def __repr__(self):
		"""Returns Analyser class full object."""
		return "Analyser({}".format(self.__dict__)

	def run(self):
		checksums = set()
		d = {}
	    # Using the Biopython fasta parse we can read our fasta input
		for record in SeqIO.parse(self.input_file, "fasta"):
			checksum = seguid(record.seq)
			if checksum in checksums:
				model_type_id = record.description.split("|")[1].split(":")[-1].strip()
				if model_type_id in ["40292"]:
					logger.warning("Duplicates {} -> {} ~ {} - [{}]".format(record.id, record.description.split("|")[-1], d[checksum].split("|")[-1], model_type_id))
				continue
			d[checksum] = record.description
			checksums.add(checksum)

