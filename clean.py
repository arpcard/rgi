import os
import sys
import filepaths
import argparse

script_path = filepaths.determine_path()
working_directory = os.getcwd()
path = script_path+"/"

app_files = [".gitignore","_docs","_tmp","_db","mgm"]

# clean other files left over
def _clean():
	import glob
	files = glob.glob(path+"*")
	for f in files:
		if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["adraft","xml","fsa","draft","pyc","log"]:
			os.remove(f)

		if os.path.isdir(f) == False:
			if os.path.isfile(f) == True and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["py"]:
				pass
			else:
				if os.path.isfile(f):
					print>>sys.stderr, "Remove: " + str(f)
					os.remove(f)

    # clean db files
	db_files = glob.glob(path+"_db/*")
	for dbfile in db_files:
		if os.path.isfile(dbfile):
			os.remove(dbfile)
	print>>sys.stderr, "Cleaned directory: "+path+"_db"

    # clean tmp files
	tmp_files = glob.glob(path+"_tmp/*")
	for tempfile in tmp_files:
		if os.path.isfile(tempfile):
			os.remove(tempfile)
	print>>sys.stderr, "Cleaned directory: "+path+"_tmp"


#remove temporary file
def main():
	if os.path.isfile(path+"contigToPro.fasta"):
		os.remove(path+"contigToPro.fasta")
	if os.path.isfile(path+"read.fsa"):
		os.remove(path+"read.fsa")
	if os.path.isfile(path+"contig.fsa"):
		os.remove(path+"contig.fsa")
	if os.path.isfile(path+"contigToORF.fsa"):
		os.remove(path+"contigToORF.fsa")
	if os.path.isfile(path+"blastRes.xml"):
		os.remove(path+"blastRes.xml")
	if os.path.isfile(path+"encodedHeaders.json"):
		os.remove(path+"encodedHeaders.json")
	if os.path.isfile(path+"predictedGenes.json"):
		os.remove(path+"predictedGenes.json")
	if os.path.isfile(path+"blastpjson"):
		os.remove(path+"blastpjson")
	if os.path.isfile(path+"proteindb.fsa"):
		os.remove(path+"proteindb.fsa")
	if os.path.isfile(path+"dnadb.fsa"):
		os.remove(path+"dnadb.fsa")
	if os.path.isfile(path+"protein.db.phr"):
		os.remove(path+"protein.db.phr")
		os.remove(path+"protein.db.pin")
		os.remove(path+"protein.db.psq")
	if os.path.isfile(path+"protein.db.dmnd"):
		os.remove(path+"protein.db.dmnd")
	if os.path.isfile(path+"temp.fsa"):
		os.remove(path+"temp.fsa")
	if os.path.isfile(path+"dna.db.nhr"):
		os.remove(path+"dna.db.nhr")
		os.remove(path+"dna.db.nin")
		os.remove(path+"dna.db.nsq")
	if os.path.isfile(path+"Report.json"):
		os.remove(path+"Report.json")
	if os.path.isfile(path+"ReportFormatted.json"):
		os.remove(path+"ReportFormatted.json")
	if os.path.isfile(path+"draft"):
		os.remove(path+"draft")
	if os.path.isfile(path+"adraft"):
		os.remove(path+"adraft")
	if os.path.isfile(path+"dataSummary.txt"):
		os.remove(path+"dataSummary.txt")

	_clean()

	print>>sys.stderr, "Cleaned directory: "+path

def run():
	parser = argparse.ArgumentParser(description='Removes BLAST databases created using card.json')
	args = parser.parse_args()
	main()

if __name__ == '__main__':
	run()
