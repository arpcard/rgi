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
		    return json.load(f)
    except ValueError as e:
        print('invalid json: %s' % e)
        return None # or: raise	

def main(filepath):
	if os.path.exists(filepath):
		if validateFile(filepath):
			dst = script_path+"/card.json"
			# copy new card.json file
			shutil.copyfile(filepath, dst)
		else:
			print "error reading json file"
	else:
		print "error uploading file"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Load card database json file',usage='runner load.py afile')
	parser.add_argument('afile',help='must be a card database json file')	
	args = parser.parse_args()
	main(sys.argv[1])