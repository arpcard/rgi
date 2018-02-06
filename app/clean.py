import argparse
import glob
from app.settings import *

# keep this files or directories
app_files = [".gitignore","_docs","_tmp","_db","mgm"]

# clean other files left over
def clean():
	files = glob.glob(path+"*")
	for f in files:
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
	data_files = glob.glob(data_path+"*")
	for datafile in data_files:
		if os.path.isfile(datafile) and os.path.basename(datafile) not in ["card.json", ".gitignore"]:
			logger.info("Remove: {}".format(datafile))
			os.remove(datafile)
	logger.info("Cleaned directory: {}".format(data_path))

    # clean db files
	db_files = glob.glob(path+"*")
	for dbfile in db_files:
		if os.path.isfile(dbfile) and os.path.basename(dbfile) not in [".gitignore"]:
			logger.info("Remove: {}".format(dbfile))
			os.remove(dbfile)
	logger.info("Cleaned directory: {}".format(path))

    # clean tmp files
	tmp_files = glob.glob(tmp+"*")
	for tempfile in tmp_files:
		if os.path.isfile(tempfile) and os.path.basename(tempfile) not in [".gitignore"]:
			logger.info("Remove: {}".format(tempfile))
			os.remove(tempfile)
	logger.info("Cleaned directory: {}".format(tmp))

    # clean log files
	log_files = glob.glob(logs+"*")
	for logfile in log_files:
		if os.path.isfile(logfile) and os.path.basename(logfile) not in [".gitignore"]:
			logger.info("Remove: {}".format(logfile))
			os.remove(logfile)
	logger.info("Cleaned directory: {}".format(logs))

def clean_local():
	if os.path.exists(LOCAL_DATABASE):
		print("clean: ", LOCAL_DATABASE)
		files = glob.glob(LOCAL_DATABASE+"*")
		for f in files:
			if os.path.isfile(f) and os.path.basename(f) not in ["card.json"]:
				logger.info("Remove: {}".format(f))
				# os.remove(f)
				print("Remove: {}".format(f))
	else:
		print("Info: Local database not found at {}, nothing to clean.".format(LOCAL_DATABASE))

#remove temporary file
def main(args):
	if args.local_database == True:
		clean_local()
	else:
		clean()


	# logger.info("Cleaned directory: {}".format(path))

def create_parser():
	parser = argparse.ArgumentParser(prog="rgi clean", description="{} - {} - Clean".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
