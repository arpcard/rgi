from app.Base import BaseModel
from app.HomologModel import Homolog
from app.VariantModel import Variant
from app.OverexpressionModel import Overexpression
from app.RrnaModel import Rrna
from app.Blast import Blast
from app.Database import Database
from app.ConvertRGIJsonToTSV import ConvertJsonToTSV
from app.settings import *

import hashlib
import multiprocessing

class Filter(BaseModel):

	"""This class takes in blast xml file and card.json file and producess perfect strict paradigm for RGI """
	def __init__(self, input_type, loose, input_sequence, xml_file, card_json, input_file, output_file, num_threads ,rgi_obj=None):
		self.input_type = input_type
		self.xml_file = xml_file
		self.card_json = card_json
		self.input_filename = input_file
		self.input_sequence = input_sequence
		self.loose = loose
		self.blast_results = {}
		self.rna_results = {}
		self.rgi_obj = rgi_obj
		self.num_threads = num_threads
		self.working_directory = rgi_obj.working_directory

		if output_file == None:
			f_path, f_name = os.path.split(input_file)
			self.output_file = os.path.join(f_path,"{}.Report.json".format(f_name))
		else:
			self.output_file = output_file

	def __repr__(self):
		"""Returns Filter class full object."""
		return "Filter({}".format(self.__dict__)

	def getHashName(self, name):
		m = hashlib.md5()
		t = time.gmtime()
		m.update(name + str(t))
		return m.hexdigest()

	def worker(self, model_type):
		logger.info("{}_worker started...".format(model_type))
		# results = {}
		try:
			if model_type == "homolog":
				obj = Homolog(self.input_type, self.loose, self.input_sequence, self.xml_file, self.working_directory, self.rgi_obj.local_database, self.rgi_obj.include_nudge)
			if model_type == "variant":
				obj = Variant(self.input_type, self.loose, self.input_sequence, self.xml_file, self.working_directory, self.rgi_obj.local_database, self.rgi_obj.include_nudge)
			if model_type == "overexpression":
				obj = Overexpression(self.input_type, self.loose, self.input_sequence, self.xml_file, self.working_directory, self.rgi_obj.local_database, self.rgi_obj.include_nudge)
			results = obj.run()
			logger.info("save {} results...".format(model_type))
			file_name = os.path.basename(self.input_sequence)
			with open(os.path.join(self.working_directory,"{}.temp.{}.json".format(file_name, model_type)), 'w') as fout:
				fout.write(json.dumps(results))

		except Exception as e:
			logger.warning("Exception: {} -> {} -> model_type: {}".format(type(e), e, model_type))

	def async_func(self):
		p = multiprocessing.Process(target=self.worker, args=("homolog",))
		p2 = multiprocessing.Process(target=self.worker, args=("variant",))
		p3 = multiprocessing.Process(target=self.worker, args=("overexpression",))
		p4 = multiprocessing.Process(target=self.process_rrna, args=("rrna",))
		# logger.debug("{} -> {}".format(p.pid, p.name))
		# logger.debug("{} -> {}".format(p2.pid, p2.name))
		# logger.debug("{} -> {}".format(p3.pid, p3.name))
		# logger.debug("{} -> {}".format(p4.pid, p4.name))
		p.start()
		p2.start()
		p3.start()
		p4.start()

	def prepare_output(self):
		"""
		Read all json into one json results file
		"""
		logger.info("prepare output(s) for input: {}".format(self.input_sequence))
		file_name = os.path.basename(self.input_sequence)
		obj=ConvertJsonToTSV(self.output_file, \
			os.path.join(self.working_directory,"{}.temp.{}.json".format(file_name, "homolog")), \
			os.path.join(self.working_directory,"{}.temp.{}.json".format(file_name, "variant")), \
			os.path.join(self.working_directory,"{}.temp.{}.json".format(file_name, "overexpression")), \
			os.path.join(self.working_directory,"{}.temp.{}.json".format(file_name, "rrna"))
			)
		# combine 3 json files			
		obj.combine_jsons()	
		# write tsv
		obj.run()

	def cleanup(self):
		self.rgi_obj.clean_files()

	def process_xml_file(self):
		""" This function is used to process blast xml file """
		model_thread = multiprocessing.Process(target=self.async_func, args=())
		model_thread.start()
		model_thread.join()
		prepare_output_thread = multiprocessing.Process(target=self.prepare_output, args=())
		prepare_output_thread.start()
		prepare_output_thread.join()
		cleanup_thread = multiprocessing.Process(target=self.cleanup, args=())
		cleanup_thread.start()

	def process_rrna(self, model_type="rrna"):
		logger.info("rRNA process: {}".format(self.input_type))
		if self.input_type == "protein":
			logger.info("Skip rRNA...")
		else:
			logger.info("rRNA process")
			self.format_fasta()
			""" Cleans rRNA model previous result and temporal files"""
			self.file_name = os.path.basename(self.input_sequence)
			d, x = self.create_db_query()
			rrna_obj = Rrna(self.input_sequence, self.output_file, d, x, self.loose, self.rgi_obj.local_database, self.rgi_obj.include_nudge)
			res = rrna_obj.run()


			file_name = os.path.basename(self.input_sequence)
			with open(os.path.join(self.working_directory,"{}.{}.json".format(file_name, model_type)), 'w') as fout:
				fout.write(json.dumps(res))			

			# with open(os.path.splitext(self.output_file)[0]+".{}.json".format(model_type), 'w') as fout:
			# 	fout.write(json.dumps(res))

			logger.info("rRNA process Done.")

	def create_db_query(self):
		logger.info("create_db_query")
		# make_custom_db(self, in_file, out_file, db_type="diamond")
		in_file = self.input_sequence
		f_path, f_name = os.path.split(self.input_sequence)
		out_file = os.path.join(self.working_directory, "{}.db".format(f_name))
		xml_file = os.path.join(self.working_directory,"{}.blastRes.rrna.xml".format(f_name))
		logger.info("DB from user query")
		db_obj = Database(self.rgi_obj.local_database)
		db_obj.make_custom_db(in_file, out_file)
		self.blast_reference_to_db_query(out_file, xml_file)
		return out_file, xml_file

	def blast_reference_to_db_query(self, db, xml_file):
		logger.info("blast_reference_to_db_query")
		# blast all rrna db against query db
		rrna_db_fasta = os.path.join(self.rgi_obj.db, "rnadb.fsa")
		blast_obj = Blast(rrna_db_fasta, program='blastn', output_file=xml_file, local_database=self.rgi_obj.local_database, num_threads=self.num_threads)
		blast_obj.run_custom(db)

	def format_fasta(self):
		f_path, f_name = os.path.split(self.input_sequence)
		temp_file = os.path.join(self.working_directory, "{}.temp".format(f_name))
		with open(temp_file, 'w') as fout:
			for record in SeqIO.parse(self.input_sequence, 'fasta'):
				header = record.id
				seq = record.seq
				fout.write(">{}\n{}\n".format(header, seq))
		self.input_sequence = temp_file

	def encode_header(self,name):
		return hashlib.md5(name.encode('utf-8')).hexdigest()

	def write_output(self):
		file_name = os.path.basename(self.input_sequence)
		logger.info(self.output_file)
		with open(self.output_file, 'w') as rrna_js:
			rrna_js.write(json.dumps(self.rna_results))

	def run(self):
		if(os.path.exists(self.xml_file)):
			self.process_xml_file()
		else:
			logger.error("missing blast xml file({}). Please check if input_type: '{}' correspond with input file: '{}' or use '--low_quality' flag for short contigs to predicts partial genes." \
					.format(self.xml_file, self.input_type, self.input_sequence))
