import os, sys
import shutil
import filepaths
import json
import argparse

# this script is used to load new card.json file to system wide package
script_path = filepaths.determine_path()

def validateFile(filename):
    try:
		with open(filename) as f:
			out = json.load(f)
		f.close()
		return out
		
    except ValueError as e:
        print>>sys.stderr, ('[error] invalid json: %s' % e)
        return None # or: raise	

def main(args):
	filepath = args.afile
	if os.path.exists(filepath):
		if validateFile(filepath):
			dst = script_path+"/_data/card.json"
			try:
				# copy new card.json file
				shutil.copyfile(filepath, dst)
				print>>sys.stderr, "[success] file copied ok"
			except Exception as e:
				print e
				print>>sys.stderr, "[error] failed to copy json file"
		else:
			print>>sys.stderr, "[error] failed to read json file"
	else:
		print>>sys.stderr,"[error] failed to upload file"

def run():
	parser = argparse.ArgumentParser(description='Load card database json file')
	parser.add_argument('-i', '--afile',help='must be a card database json file')
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()
