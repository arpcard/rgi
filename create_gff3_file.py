import json
import csv
import sys

'''
Column 1: "seqid"
-----------------
	The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped whitespace and must not begin with an unescaped ">".

Column 2: "source"
-----------------
	The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type that is a subclass of the type in the type column.

Column 3: "type"
-----------------
	The type of the feature (previously called the "method"). This is constrained to be either: (a)a term from the "lite" version of the Sequence Ontology - SOFA, a term from the full Sequence Ontology - it must be an is_a child of sequence_feature (SO:0000110) or (c) a SOFA or SO accession number. The latter alternative is distinguished using the syntax SO:000000.

Columns 4 & 5: "start" and "end"
--------------------------------
	The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be less than or equal to end is satisfied by making end = the position of the end + the length of the landmark feature.

	For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark.

Column 6: "score"
-----------------
	The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features.

Column 7: "strand"
-----------------
	The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.

Column 8: "phase"
-----------------
	For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. In other words, 
		a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, 
		a phase of "1" indicates that the next codon begins at the second base of this region, 
		and a phase of "2" indicates that the codon begins at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3.

	For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field.

	The phase is REQUIRED for all CDS features.

Column 9: "attributes"
----------------------
	A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.

	These tags have predefined meanings:

	ID
		Indicates the ID of the feature. IDs for each feature must be unique within the scope of the GFF file. In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID collectively represent a single feature.
	Name
		Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
	Alias
		A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.
	Parent
		Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship.
	Target
		Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
	Gap
		The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is taken from the CIGAR format described in the Exonerate documentation. See "THE GAP ATTRIBUTE" for a description of this format.
	Derives_from
		Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See "PATHOLOGICAL CASES" for further discussion.
	Note
		A free text note.
	Dbxref
		A database cross reference. See the section "Ontology Associations and Db Cross References" for details on the format.
	Ontology_term
		A cross reference to an ontology term. See the section "Ontology Associations and Db Cross References" for details.
	Is_circular
		A flag to indicate whether a feature is circular. See extended discussion below.
'''

def create_gff3(afile,orf):
	try:
		with open(afile, 'r') as f:
			data = json.load(f)
		f.close()
	
	except ValueError:
		print>>sys.stderr, "convertJsonToTSV expects a file contains a VALID JSON string."
		exit()
	
	with open("gff3.txt", "w") as af:
		writer = csv.writer(af, delimiter='\t', dialect='excel')
		writer.writerow(['##gff-version 3'])
		#writer.writerow(["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
		for item in data:
			minevalue = False
			maxpercent = 0
			startCompare = False
			minARO = 0
			minARO_id = 0
			hitID = ""
			geneID = ""
			AROlist = {}

			if item not in ["_metadata","data_type"]:
				geneID = item
				for it in data[item]:

					idenPercent = float(data[item][it]["max-identities"]) / len(data[item][it]["query"])
					# sort results by minimum e-value and maximum percent identity
					AROlist[str(data[item][it]["ARO_name"])] = data[item][it]["ARO_accession"]

					if startCompare:
						if minevalue > data[item][it]["evalue"] and float(maxpercent) < float(idenPercent):
							minevalue = data[item][it]["evalue"]
							maxpercent = float(idenPercent)
							minARO = data[item][it]["ARO_name"]
							minARO_id = data[item][it]["ARO_accession"]
							SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]

							if "orf_prot_sequence" in data[item][it]:
								predictedProtein = data[item][it]["orf_prot_sequence"]
							if "hsp_num:" in it:
								hitID = it

					else:
						startCompare = True
						minevalue = data[item][it]["evalue"]
						maxpercent = float(idenPercent)
						minARO = data[item][it]["ARO_name"]
						minARO_id = data[item][it]["ARO_accession"]
						SequenceFromBroadStreet = data[item][it]["SequenceFromBroadStreet"]

						if "orf_prot_sequence" in data[item][it]:
							predictedProtein = data[item][it]["orf_prot_sequence"]
						if "hsp_num:" in it:
							hitID = it

		 		_seqid = str(geneID) + "|" + str(hitID)
				_source = "RGI:3.1.2"	
				_type = "CDS"
				_phase = "."
				_attributes = "name=" + str(minARO)+";alias=" + "ARO:"+str(minARO_id)
				_score = str(minevalue)

				if orf == "genemark":
					_start = str(int(findnthbar(item, 4))-1)
					_end = str(int(findnthbar(item, 5))-1)
					_strand = str(findnthbar(item, 3))				
				else:
					_start = str(int(findnthbar2(item, 1))-1)
					_end = str(int(findnthbar2(item, 2))-1)
					_strand = str(findnthbar2(item, 3))

				writer.writerow([_seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes])
	af.close()

#output the information particular field from alignment.Title by splicing it by '|'
def findnthbar(bunchstr, start):

	barc = 0
	over = start+1
	temp = ""

	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if barc == start:
			if eachc == '|':
				pass
			else:
				temp += eachc
		if barc == over:
			break		

	return temp

#output the information particular field from alignment.Title by splicing it by '#'
def findnthbar2(bunchstr, n):
	arr = bunchstr.split("#")
	if n < len(arr):
		# gene id
		if n == 1 and arr[n]:
			return int(arr[n])
		elif n == 2:
			return int(arr[n])
		elif n == 3:
			 if int(arr[n]) == 1:
			 	# positive
			 	return "+"
			 else:
			 	# neg
			 	return "-"
		else:
			return arr[n]
	else:
		return ""

def main(afile,orf):
	create_gff3(afile,orf)


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])
	"""required: sys.argv[1] must be a json file"""