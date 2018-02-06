import shutil
import argparse
from app.settings import *

# this script is used to load new card.json file to system wide package or local

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
			# path to save card.json file and database files
			if args.local_database == True:
				db = LOCAL_DATABASE
				# create directory if it doesn't exist
				if not os.path.exists(LOCAL_DATABASE):
					os.makedirs(LOCAL_DATABASE)
			else:
				db = data_path
			try:
				# copy new card.json file
				shutil.copyfile(filepath, os.path.join(db, "card.json"))
				logger.info("[success] file copied ok")
				print("[success] file copied ok")
			except Exception as e:
				logger.warning("failed to copy json file: {}".format(e))
		else:
			logger.warning("[error] failed to read json file")
	else:
		logger.warning("[error] failed to upload file")

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('-i', '--afile',help='must be a card database json file')
	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()
