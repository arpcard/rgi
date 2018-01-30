import shutil
import argparse
from app.settings import *

# this script is used to load new card.json file to system wide package

def validateFile(filename):
	try:
		with open(filename) as f:
			out = json.load(f)
			return out
	except ValueError as e:
		logger.error("Invalid json: {}".format(e))
		return None # or: raise

def main(args):
	filepath = args.afile
	if os.path.exists(filepath):
		if validateFile(filepath):
			try:
				# copy new card.json file
				shutil.copyfile(filepath, os.path.join(data_path, "card.json"))
				logger.info("[success] file copied ok")
				print("[success] file copied ok")
			except Exception as e:
				logger.warning("failed to copy json file: {}".format(e))
		else:
			logger.warning("[error] failed to read json file")
	else:
		logger.warning("[error] failed to upload file")

def run():
	parser = argparse.ArgumentParser(description='Load card database json file')
	parser.add_argument('-i', '--afile',help='must be a card database json file')
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()
