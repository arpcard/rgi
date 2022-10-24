import argparse
import glob
import shutil
from app.settings import *

# clean other files left over
def clean():
	files = glob.glob(os.path.join(path,"*"))
	for f in files:
		remove_directory(f)
		if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["adraft","xml","fsa","draft","pyc","log"]:
			os.remove(f)

		if os.path.isdir(f) == False:
			if os.path.isfile(f) == True and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["py","md"]:
				pass
			else:
				if os.path.isfile(f):
					logger.info("Remove: {}".format(f))
					os.remove(f)

    # clean data files
	data_files = glob.glob(os.path.join(data_path,"*"))
	for datafile in data_files:
		if os.path.isfile(datafile) and os.path.basename(datafile) not in ["card.json", ".gitignore","__init__.py"]:
			logger.info("Remove: {}".format(datafile))
			os.remove(datafile)
	logger.info("Cleaned directory: {}".format(data_path))

    # clean db files
	db_files = glob.glob(os.path.join(path,"*"))
	for dbfile in db_files:
		if os.path.isfile(dbfile) and os.path.basename(dbfile) not in [".gitignore","__init__.py"]:
			logger.info("Remove: {}".format(dbfile))
			os.remove(dbfile)
	logger.info("Cleaned directory: {}".format(path))

def clean_local():
	if os.path.exists(LOCAL_DATABASE):
		logger.info("clean: {}".format(LOCAL_DATABASE))
		files = glob.glob(os.path.join(LOCAL_DATABASE,"*"))
		for f in files:
			if os.path.isfile(f) and os.path.basename(f) not in ["card.json"]:
				logger.info("Remove: {}".format(f))
				os.remove(f)
			# remove bwt directory
			remove_directory(f)
	else:
		logger.warning("Local database not found at {}, nothing to clean.".format(LOCAL_DATABASE))

def remove_directory(directory_path):
	if os.path.basename(directory_path) in ["bwt"] and os.path.isdir(directory_path):
		try:
			logger.info("Remove dir: {}".format(directory_path))
			shutil.rmtree(directory_path)
		except OSError as e:
			logger.error("{} - {}" % (e.filename, e.strerror))

#remove temporary file
def main(args):
	if args.debug:
		logger.setLevel(10)
	if args.local_database == True:
		clean_local()
	else:
		clean()

	# logger.info("Cleaned directory: {}".format(path))

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi clean", description="{} - {} - Clean".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
