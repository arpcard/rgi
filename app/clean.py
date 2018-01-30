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

#remove temporary file
def main():
	clean()
	# logger.info("Cleaned directory: {}".format(path))

def run():
	parser = argparse.ArgumentParser(description='Removes BLAST databases created using card.json')
	args = parser.parse_args()
	main()

if __name__ == '__main__':
	run()
