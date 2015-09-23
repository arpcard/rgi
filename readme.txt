Documentation of rgi.py

Before you run RGI script, make sure you have installed all other bioinformatics tool already:

* MetaGeneMark http://exon.gatech.edu/GeneMark/license_download.cgi
* BLAST ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* Biopython http://biopython.org/DIST/docs/install/Installation.html#sec12
* Download the database - card.json from Downloads in our website

Open a terminal, type: 

> python rgi.py inputSequenceType inputSequence 

eg. python rgi.py protein query.fasta

Currently, inputSequenceType could be one of 'contig', 'protein' or 'read'.
inputSequence is a file in a specific text format in bioinformatics. 

1. 'Contig' means that inputSequence is a DNA sequence stored in a FastA file. It can be preprocessed to find an open reading frame (orf) in a sequence or a complete genome sequence.
2. 'Protein', as its name suggests, requires a Fasta file with protein sequence. 
3. 'Read', however, is in FastQ format, which is a completely raw DNA data from experiments. 

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
 