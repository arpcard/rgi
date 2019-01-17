from app.settings import *
import glob, subprocess

class Baits(object):
	"""Class to create Baits object and get melting points for DNA probes."""
	def __init__(self,input_file, output_file, filter_temperature, clean, debug):
		"""Creates Baits object"""
		self.input_file = input_file
		self.output_file = output_file
		self.filter_temperature = filter_temperature

		o_f_path, o_f_name = os.path.split(os.path.abspath(self.output_file))

		self.all_output = os.path.join(o_f_path, self.output_file + ".json")
		self.filtered_output = os.path.join(o_f_path, self.output_file + "_filtered_{}.json".format(self.filter_temperature))
		self.clean = clean
		self.debug = debug
		
		if self.debug:
			logger.setLevel(10)

	def __repr__(self):
		"""Returns Baits class full object."""
		return "Baits({}".format(self.__dict__)

	def check_melt_program(self):
		logger.info("check for melt.pl ...")
		output = ""

		try:
			output = subprocess.check_output("which melt.pl", shell=True)
		except Exception as e:
			output = e

		if "returned non-zero exit status 1" in str(output):
			logger.error("Missing program melt.pl. Please install melt.pl from: http://unafold.rna.albany.edu/?q=unafold-man-pages")
			exit()

	def run(self):

		self.check_melt_program()
		logger.info("calculating baits melting temperature and entropy ...")
		os.system('{program} -n RNA {input} > {output_file}' \
					.format(
						program="melt.pl", 
						input=self.input_file,
						output_file=self.output_file
					)
				)

		# parse output file
		probes = {}
		temps = []

		with open(self.output_file, 'r') as in_file:
			data = in_file.readlines()
			startrecord = False
			count = 1
			for eachline in data:
				if eachline == '':
					pass
				elif eachline[0:16] == "Calculating for ":
					probes[count] = {
						"id": eachline.split(", t = ")[0].split("Calculating for ")[-1],
						"raw": eachline,  
						"change_in_gibbs_free_energy (dG)": 0.0, 
						"change_in_enthalpy (dH)": 0.0, 
						"change_in_entropy (dS)": 0.0,
						"melting_temperature (Tm)": 0.0,
					}
					count = count + 1

			c = 1

			for line in data:
				if line == '':
					pass
				elif line[0:11] == "dG	dH	dS	Tm":
					startrecord = True
				elif startrecord:
					temp = line.strip("\n").split("\t")
					probes[c]["change_in_gibbs_free_energy (dG)"] = float(temp[0])
					probes[c]["change_in_enthalpy (dH)"] = float(temp[1])
					probes[c]["change_in_entropy (dS)"] = float(temp[2])
					probes[c]["melting_temperature (Tm)"] = float(temp[3])
					c = c + 1

			logger.info("saved results in {}".format(self.all_output))
			with open(self.all_output, 'w') as out_file:
				json.dump(probes, out_file)

		# filter by desired melting temperature
		self.melt_filter()

		# remove temporary files
		if self.clean == True:
			o_f_path, o_f_name = os.path.split(os.path.abspath(self.output_file))
			files = glob.glob(os.path.join(o_f_path,"*"))
			for f in files:
				if os.path.isfile(f) == True and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["run","plot","ext","ct","dG"]:
					logger.info("Remove: {}".format(f))
					os.remove(f)


	def melt_filter(self):
		filtered_results = {}
		j = 1
		with open(self.all_output, 'r') as jfile:
			data = json.load(jfile)
			for i in data:
				if data[i]["melting_temperature (Tm)"] >= float(self.filter_temperature):
					filtered_results[j] = data[i]
					j = j + 1

		logger.info("saved filtered results (tempature by {}) in {}".format(self.filter_temperature, self.filtered_output))

		with open(self.filtered_output, 'w') as out_file:
			json.dump(filtered_results, out_file)




