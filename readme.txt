Documentation of rgi.py

Before you run the RGI scripts, make sure you have installed needed external tools:

* MetaGeneMark http://exon.gatech.edu/GeneMark/license_download.cgi

- Tested with MetaGeneMark v3.26 on linux 64 and Mac OS X
- After downloading MetaGeneMark copy mgm directory into release-rgi folder
- Follow the INSTALL instructions inside the mgm
- Change directory to release-rgi
- Test MetaGeneMark using the following command:

> ./mgm/gmhmmp (There should be list of commands)

* BLAST ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

- Tested with BLAST 2.2.28 and BLAST 2.2.31+ on linux 64 and Mac OS X

* You can alson run the following command to install blast. This will only install version 2.2.28

> sudo apt-get install ncbi-blast+

- Test blast install with the following command:

> makeblastdb

* Biopython http://biopython.org/DIST/docs/install/Installation.html#sec12
* Run the following command to install Bio-python

> sudo apt-get install python-biopython

* Download the database - card.json from Downloads on the CARD website (a copy may be included with this release)

Running RGI:

Open a terminal, type: 

> python rgi.py inputSequenceType inputSequence 

eg. python rgi.py protein query.fasta

Currently, inputSequenceType could be one of 'contig', 'protein' or 'read'.

1. 'contig' means that inputSequence is a DNA sequence stored in a FASTA file, presumably a complete genome or assembly contigs. RGI will predict ORFs de novo and predict resistome using a combination of BLASTP against the CARD data, curated cut-offs, and SNP screening.

2. 'protein', as its name suggests, requires a FASTA file with protein sequences. As above, RGI predict resistome using a combination of BLASTP against the CARD data, curated cut-offs, and SNP screening.

3. 'read' expects raw FASTQ format nucleotide data and predicts resistome using a combination of BLASTX against the CARD data, curated cut-offs, and SNP screening. This is an experimental tool and we have yet to adjust the CARD cut-offs for BLASTX.  We will be exploring other metagenomics or FASTQ screening methods. Note that RGI does not perform any pre-processing of the FASTQ data (linker trimming, etc).

RGI Output will produce a detailed JSON file: Report.json

The JSON is as follows (example shows only one hit):

- gene_71|gi|378406451|gb|JN420336.1| Klebsiella pneumoniae plasmid pNDM-MAR, complete sequence: {
	// Hit 1
	gnl|BL_ORD_ID|39|hsp_num:0: {
		SequenceFromBroadStreet: "MRYIRLCIISLLATLPLAVHASPQPLEQIKQSESQLSGRVGMIEMDLASGRTLTAWRADERFPMMSTFKVVLCGAVLARVDAGDEQLERKIHYRQQDLVDYSPVSEKHLADGMTVGELCAAAITMSDNSAANLLLATVGGPAGLTAFLRQIGDNVTRLDRWETELNEALPGDARDTTTPASMAATLRKLLTSQRLSARSQRQLLQWMVDDRVAGPLIRSVLPAGWFIADKTGASKRGARGIVALLGPNNKAERIVVIYLRDTPASMAERNQQIAGIGAA",
		"orf_start": 67822,
		"ARO_name": "SHV-12",
		"type_match": "Loose",
		"query": "INDWRLDYNECRPHSSLNYLTPAEFAAGWRN",
		"evalue": 3.82304,
		"max-identities": 10,
		"orf_strand": "-",
		"bit-score": 24.6386,
		"cvterm_id": "35914",
		"sequenceFromDB": "LDRWETELNEALPGDARDTTTPASMAATLRK",
		"match": "++ W  + NE  P  + +  TPA  AA  R ",
		"model_id": "103",
		"orf_From": "gi|378406451|gb|JN420336.1| Klebsiella pneumoniae plasmid pNDM-MAR, complete sequence",
		"pass_evalue": 1e-100,
		"query_end": 68607,
		"ARO_category": {
		    "36696": {
			    "category_aro_name": "antibiotic inactivation enzyme",
			    "category_aro_cvterm_id": "36696",
			    "category_aro_accession": "3000557",
			    "category_aro_description": "Enzyme that catalyzes the inactivation of an antibiotic resulting in resistance.  Inactivation includes chemical modification, destruction, etc."
			},
			"36268": {
		        "category_aro_name": "beta-lactam resistance gene",
		        "category_aro_cvterm_id": "36268",
		        "category_aro_accession": "3000129",
		        "category_aro_description": "Genes conferring resistance to beta-lactams."
		    }
		},
		"ARO_accession": "3001071",
		"query_start": 68515,
		"model_name": "SHV-12",
		"model_type": "model-blastP",
		"orf_end": 68646
	},
	...
	// Hit 2
	...
	// Hit 3
	...
}

Run the following command to get the Tab Delimited output

> python convertJsonToTSV.py ./Report.json

This outputs a tab-delimited text file: dataSummary.txt

The tab-output is as follows:

ORF_ID	CONTIG	START	STOP	ORIENTATION	CUT_OFF	Best_Hit_evalue	Best_Hit_ARO	Best_Identites	ARO	ARO_name	Model_type	SNP	AR0_category	bit_score

Example:

ORF_ID: gene_71

CONTIG: gi|378406451|gb|JN420336.1| Klebsiella pneumoniae plasmid pNDM-MAR, complete sequence

START: 67822

STOP: 68646

ORIENTATION: -

CUT_OFF: Loose

Best_Hit_evalue: 0.183312

Best_Hit_ARO: vanWB

Best_Identites: 0.35714285714285715

ARO: ARO:3001071, ARO:3001136, ARO:3002422, ARO:3001202, ARO:3001336, ARO:3001075, ARO:3001108, ARO:3001122, ARO:3001340, ARO:3001083, ARO:3001199, ARO:3001182, ARO:3001124, ARO:3001073, ARO:3001147, ARO:3001099, ARO:3001144, ARO:3001088, ARO:3001145, ARO:3001155, ARO:3001198, ARO:3001183, ARO:3001148, ARO:3001176, ARO:3001068, ARO:3001134, ARO:3001175, ARO:3001112, ARO:3002909, ARO:3001121, ARO:3001131, ARO:3001159, ARO:3001177, ARO:3001060, ARO:3001192, ARO:3002472, ARO:3001344, ARO:3001116, ARO:3001345, ARO:3001178, ARO:3001089, ARO:3001191, ARO:3001125, ARO:3001173, ARO:3001197, ARO:3001119, ARO:3001117, ARO:3001157, ARO:3001168, ARO:3001086, ARO:3001337, ARO:3001141, ARO:3001352, ARO:3001087, ARO:3001137, ARO:3001082, ARO:3001167, ARO:3001158, ARO:3001187, ARO:3001204, ARO:3001132, ARO:3002424, ARO:3001152, ARO:3001149, ARO:3001070, ARO:3001194, ARO:3001185, ARO:3001109, ARO:3001067, ARO:3001362, ARO:3003156, ARO:3001129, ARO:3001139, ARO:3001093, ARO:3001189, ARO:3001101, ARO:3001200, ARO:3001079, ARO:3001065, ARO:3001085, ARO:3001338, ARO:3001100, ARO:3001146, ARO:3001184, ARO:3001127, ARO:3001105, ARO:3001181, ARO:3001356, ARO:3001196, ARO:3001113, ARO:3001074, ARO:3001361, ARO:3001072, ARO:3003152, ARO:3001091, ARO:3001169, ARO:3001172, ARO:3001195, ARO:3001186, ARO:3001203, ARO:3001126, ARO:3001114, ARO:3001150, ARO:3001151, ARO:3001171, ARO:3001156, ARO:3001111, ARO:3001110, ARO:3002964, ARO:3002475, ARO:3001106, ARO:3001120, ARO:3001061, ARO:3001092, ARO:3001103, ARO:3001140, ARO:3001064, ARO:3001095, ARO:3001098, ARO:3001174, ARO:3001364, ARO:3001104, ARO:3001188, ARO:3001170, ARO:3001128, ARO:3001138, ARO:3001059, ARO:3001357, ARO:3001066, ARO:3001347, ARO:3001090, ARO:3001201, ARO:3001190, ARO:3001123, ARO:3001135, ARO:3001153, ARO:3001193, ARO:3001102, ARO:3001118, ARO:3001351, ARO:3001096, ARO:3001094, ARO:3001077, ARO:3001076, ARO:3001080, ARO:3001154, ARO:3001130, ARO:3001078, ARO:3003154, ARO:3001161

ARO_name: SHV-75, SHV-74, SHV-77, SHV-76, SHV-71, SHV-70, SHV-73, SHV-72, SHV-78, SHV-144, SHV-145, SHV-147, SHV-140, SHV-141, SHV-142, SHV-143, SHV-148, SHV-149, SHV-189, SHV-40, SHV-41, SHV-42, SHV-43, SHV-44, SHV-45, SHV-46, SHV-48, SHV-49, SHV-5, SHV-173, SHV-172, SHV-179, SHV-178, SHV-59, SHV-57, SHV-56, SHV-55, SHV-53, SHV-52, SHV-51, SHV-168, SHV-167, SHV-164, SHV-165, SHV-162, SHV-163, SHV-160, SHV-161, SHV-2A, SHV-28, SHV-29, SHV-22, SHV-20, SHV-21, SHV-27, SHV-24, SHV-25, vanWB, SHV-119, SHV-112, SHV-110, SHV-38, SHV-182, SHV-183, SHV-185, SHV-187, SHV-31, SHV-30, SHV-33, SHV-32, SHV-35, SHV-34, SHV-37, SHV-36, SHV-108, SHV-109, SHV-127, SHV-100, SHV-101, SHV-102, SHV-103, SHV-104, SHV-105, SHV-106, SHV-107, SHV-134, LEN-3, SHV-99, LEN-4, SHV-2, SHV-137, SHV-154, SHV-7, SHV-6, SHV-133, SHV-158, SHV-9, SHV-8, vanG, OKP-A-5, SHV-125, SHV-122, SHV-123, SHV-120, SHV-121, SHV-126, SHV-95, SHV-124, OKP-A-7, SHV-128, SHV-129, SHV-89, SHV-84, SHV-85, SHV-86, SHV-80, SHV-81, SHV-82, SHV-83, SHV-66, SHV-67, SHV-64, SHV-65, SHV-62, SHV-63, SHV-61, SHV-69, SHV-1, SHV-157, SHV-156, SHV-155, SHV-98, SHV-153, SHV-152, SHV-151, SHV-150, SHV-93, SHV-92, SHV-97, SHV-96, SHV-159, SHV-94, SHV-13, SHV-12, SHV-11, SHV-16, SHV-15, SHV-14, SHV-19, SHV-18

Model_type: model-blastP

SNP: n/a

AR0_category: antibiotic inactivation enzyme, antibiotic resistance gene cluster, cassette, or operon, beta-lactam resistance gene, gene conferring antibiotic resistance via molecular bypass, glycopeptide resistance gene

bit_score: 24.6386, 24.6386, 23.483, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 26.9498, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 23.483, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 23.8682, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 23.483, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 26.1794, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.2534, 23.483, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 23.8682, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 28.8758, 23.483, 24.6386, 24.6386, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 24.6386, 23.8682, 24.6386, 24.6386, 24.2534, 24.6386, 24.6386, 24.2534, 24.6386, 24.2534, 24.6386, 23.8682
					
Windows:

* Download and install vmware at : https://my.vmware.com/web/vmware/free#desktop_end_user_computing/vmware_workstation_player/12_0
* Download ubuntu iso at: http://www.ubuntu.com/download/desktop
* Create a virtual (Enable virtualization technology)

Beta Testing Notes:

* RGI file (rgi.py) is from build api-master-8dcc61351daa02b5c1da8a1b119f26032f668c81
* Database is created once the rgi.py is run. Core files are:

- card.json
- clean.py
- contigToProteins.py
- convertJsonToTSV.py
- fqToFsa.py
- readme.txt
- rgi.py
- ./mgm/ (This directory contains GeneMark software and license_download. Download this software separately)

* If new card.json is available. Replace card.json in this directory.
* Then run clean.py to clean directory.

> python clean.py

Contact Us:

For help please contact the following awesome people:

* Dr. Andrew McArthur <mcarthua@mcmaster.ca>
* Amos Raphenya <raphenar@mcmaster.ca>
* Pearl Guo <p7guo@uwaterloo.ca>
* Justin Jia <jiabf@mcmaster.ca>
 