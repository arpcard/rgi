import os
import sys

#remove temporary file
def main():
	
	if os.path.isfile("contigToPro.fasta"):
		os.remove("contigToPro.fasta")
	if os.path.isfile("read.fsa"):
		os.remove("read.fsa")
	if os.path.isfile("contig.fsa"):
		os.remove("contig.fsa")
	if os.path.isfile("contigToORF.fsa"):
		os.remove("contigToORF.fsa")
	if os.path.isfile("blastRes.xml"):
		os.remove("blastRes.xml")
	if os.path.isfile("blastpjson"):
		os.remove("blastpjson")
	if os.path.isfile("proteindb.fsa"):
		os.remove("proteindb.fsa")
	if os.path.isfile("dnadb.fsa"):
		os.remove("dnadb.fsa")
	if os.path.isfile("protein.db.phr"):
		os.remove("protein.db.phr")
		os.remove("protein.db.pin")
		os.remove("protein.db.psq")
	if os.path.isfile("temp.fsa"):
		os.remove("temp.fsa")
	if os.path.isfile("dna.db.nhr"):
		os.remove("dna.db.nhr")
		os.remove("dna.db.nin")
		os.remove("dna.db.nsq")

	print>>sys.stderr, "Cleaned directory"

if __name__ == '__main__':
	main()