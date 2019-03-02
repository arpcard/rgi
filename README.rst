.. image:: https://travis-ci.org/arpcard/rgi.svg?branch=master
    :target: https://travis-ci.org/arpcard/rgi

Resistance Gene Identifier (RGI) 
--------------------------------------------

This application is used to predict resistome(s) from protein or nucleotide data based on homology and SNP models. The application uses reference data from the `Comprehensive Antibiotic Resistance Database (CARD) <https://card.mcmaster.ca/>`_.

RGI analyses can be performed via the CARD website `RGI portal <https://card.mcmaster.ca/analyze/rgi>`_, via use of a `Galaxy wrapper <https://github.com/arpcard/rgi_wrapper>`_ for the `Galaxy <https://galaxyproject.org/tutorials/g101>`_ platform, or alternatively you can `Install RGI from Conda`_ or `Run RGI from Docker`_. The instructions below discuss use of RGI at the command line, following a general overview of how RGI works for genomes, genome assemblies, proteomes, and metagenomic sequencing.

Analyzing Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes
-----------------------------------------------------------------------

If DNA sequences are submitted, RGI first predicts complete open reading frames (ORFs) using `Prodigal <https://github.com/hyattpd/Prodigal>`_ (ignoring those less than 30 bp) and analyzes the predicted protein sequences. Short contigs, small plasmids, low quality assemblies, or merged metagenomic reads should be analyzed using Prodigal's algorithms for low quality/coverage assemblies (i.e. contigs <20,000 bp) and inclusion of partial gene prediction. If the low sequence quality option is selected, RGI uses Prodigal anonymous mode for open reading frame prediction, supporting calls of partial AMR genes from short or low quality contigs.

If protein sequences are submitted, RGI skips ORF prediction and uses the protein sequences directly.

The RGI currently supports CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_ (use of BLAST or `DIAMOND <https://ab.inf.uni-tuebingen.de/software/diamond>`_ bitscore cut-offs to detect functional homologs of AMR genes), `protein variant models <https://card.mcmaster.ca/ontology/40293>`_ (for accurate differentiation between susceptible intrinsic genes and intrinsic genes that have acquired mutations conferring AMR, based on CARD's curated SNP matrices), `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_ (for detection of drug resistant rRNA target sequences), and `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_ (which detect efflux subunits associated AMR, but also highlights mutations conferring over-expression when present).

The RGI analyzes genome or proteome sequences under three paradigms: Perfect, Strict, and Loose (a.k.a. Discovery). The Perfect algorithm is most often applied to clinical surveillance as it detects perfect matches to the curated reference sequences and mutations in the CARD. In contrast, the Strict algorithm detects previously unknown variants of known AMR genes, including secondary screen for key mutations, using detection models with CARD's curated similarity cut-offs to ensure the detected variant is likely a functional AMR gene. The Loose algorithm works outside of the detection model cut-offs to provide detection of new, emergent threats and more distant homologs of AMR genes, but will also catalog homologous sequences and spurious partial hits that may not have a role in AMR. Combined with phenotypic screening, the Loose algorithm allows researchers to hone in on new AMR genes.

**Note: All Loose hits of 95% identity or better are automatically listed as Strict.**

All results are organized via the `Antibiotic Resistance Ontology <https://card.mcmaster.ca/ontology/36006>`_ classification: AMR Gene Family, Drug Class, and Resistance Mechanism.

Analyzing Metagenomic Reads
--------------------------------------------

Insert text here.

Table of Contents
-------------------------------------

- `License`_
- `Citation`_
- `Support & Bug Reports`_
- `Requirements`_
- `Install Dependencies`_
- `Install RGI from Project Root`_
- `Running RGI Tests`_
- `Help Menu and Usage`_
- `Help Menus for Subcommands`_
- `Load card.json`_
- `Check Database Version`_
- `Clean Previous or Old Databases`_
- `RGI main Usage for Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes`_
- `Running RGI main with Genome or Assembly DNA Sequences`_
- `Running RGI main with Protein Sequences`_
- `Running RGI main using GNU Parallel`_
- `Generating Heat Maps of RGI main Results`_
- `Run RGI from Docker`_
- `Install RGI from Conda`_
- `Overview of Tab-Delimited Output`_

License
--------

Use or reproduction of these materials, in whole or in part, by any non-academic organization whether or not for non-commercial (including research) or commercial purposes is prohibited, except with written permission of McMaster University. Commercial uses are offered only pursuant to a written license and user fee. To obtain permission and begin the licensing process, see the `CARD website <https://card.mcmaster.ca/about>`_.

Citation
--------

Jia et al. 2017. CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 45, D566-573. [`PMID 27789705 <https://www.ncbi.nlm.nih.gov/pubmed/27789705>`_]

Support & Bug Reports
----------------------

Please log an issue on `github issue <https://github.com/arpcard/rgi/issues>`_.

You can email the CARD curators or developers directly at `card@mcmaster.ca <mailto:card@mcmaster.ca>`_, via Twitter at `@arpcard <http://www.twitter.com/arpcard>`_.

Requirements
--------------------

- `Python 3.6 <https://www.python.org/>`_
- `NCBI BLAST 2.6.0 <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
- `six 1.7.0+ <https://bitbucket.org/gutworth/six>`_
- `zlib <https://bitbucket.org/gutworth/six>`_
- `Prodigal 2.6.3 <https://github.com/hyattpd/prodigal/wiki/Installation>`_
- `DIAMOND 0.8.36 <https://ab.inf.uni-tuebingen.de/software/diamond>`_
- `Biopython 1.60+ <https://biopython.org/>`_
- `filetype 1.0.0+ <https://pypi.org/project/filetype/>`_
- `pytest 3.0.0+ <https://docs.pytest.org/en/latest/>`_
- `mock 2.0.0 <https://pypi.org/project/mock/>`_
- `pandas 0.15.0+ <https://pandas.pydata.org/>`_
- `Matplotlib 2.1.2+ <https://matplotlib.org/>`_
- `seaborn 0.8.1+ <https://matplotlib.org/>`_
- `pyfaidx 0.5.4.1+ <https://pypi.org/project/pyfaidx/>`_
- `pyahocorasick 1.1.7+ <https://pypi.org/project/pyahocorasick/>`_
- `OligoArrayAux 3.8 <http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux>`_
- `samtools 1.9 <https://github.com/samtools/samtools>`_
- `bamtools 2.5.1 <https://github.com/pezmaster31/bamtools>`_
- `bedtools 2.27.1 <https://github.com/arq5x/bedtools2>`_
- `Jellyfish 2.2.10 <https://github.com/gmarcais/Jellyfish>`_
- `Bowtie2 2.3.4.3 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
- `BWA 0.7.17 (r1188) <https://github.com/lh3/bwa>`_

Install Dependencies
--------------------

- pip3 install six
- pip3 install biopython
- pip3 install filetype
- pip3 install pytest
- pip3 install mock
- pip3 install pandas
- pip3 install matplotlib
- pip3 install seaborn
- pip3 install pyfaidx
- pip3 install pyahocorasick

Install RGI from Project Root
-----------------------------

.. code-block:: sh

   pip3 install .

or

.. code-block:: sh

   python3 setup.py build
   python3 setup.py test
   python3 setup.py install

Running RGI Tests
-------------------
.. code-block:: sh
   
   cd tests
   pytest -v -rxs

Help Menu and Usage
----------------------

The following command will bring up RGI's main help menu:

.. code-block:: sh

   rgi --help

.. code-block:: sh

      usage: rgi <command> [<args>]
                  commands are:
                  ---------------------------------------------------------------------------------------
                  Database
                  ---------------------------------------------------------------------------------------

                  load     Loads CARD database, annotations and k-mer database
                  clean    Removes BLAST databases and temporary files
                  database Information on installed card database
                  galaxy   Galaxy project wrapper

                  ---------------------------------------------------------------------------------------
                  Genomic
                  ---------------------------------------------------------------------------------------

                  main     Runs rgi application
                  tab      Creates a Tab-delimited from rgi results
                  parser   Creates categorical JSON files RGI wheel visualization
                  heatmap  Heatmap for multiple analysis

                  ---------------------------------------------------------------------------------------
                  Metagenomic
                  ---------------------------------------------------------------------------------------
                  bwt                   Align reads to CARD and in silico predicted allelic variants
                  
                  ---------------------------------------------------------------------------------------
                  Baits validation
                  ---------------------------------------------------------------------------------------
                  tm                    Baits Melting Temperature

                  ---------------------------------------------------------------------------------------
                  Annotations
                  ---------------------------------------------------------------------------------------
                  card_annotation       Create fasta files with annotations from card.json
                  wildcard_annotation   Create fasta files with annotations from variants
                  baits_annotation      Create fasta files with annotations from baits (Experimental)
                  remove_duplicates     Removes duplicate sequences (Experimental)

                  ---------------------------------------------------------------------------------------
                  Pathogen of origin
                  ---------------------------------------------------------------------------------------
                  
                  kmer_build            Build AMR specific k-mers database used for pathogen of origin
                  kmer_query            Query sequences against AMR k-mers database to predict pathogen of origin

   Resistance Gene Identifier - <version_number>

   positional arguments:
   command     Subcommand to run

   optional arguments:
   -h, --help  show this help message and exit

   Use the Resistance Gene Identifier to predict resistome(s) from protein or
   nucleotide data based on homology and SNP models. Check
   https://card.mcmaster.ca/download for software and data updates. Receive email
   notification of monthly CARD updates via the CARD Mailing List
   (https://mailman.mcmaster.ca/mailman/listinfo/card-l)

Help Menus for Subcommands
----------------------------

Help screens for subcommands can be accessed using the -h argument, e.g.

.. code-block:: sh

      rgi load -h

Load card.json 
-------------------

To start analyses, first acquire the latest AMR reference data from CARD at `https://card.mcmaster.ca/latest/data <https://card.mcmaster.ca/latest/data>`_. CARD data can be installed at the system level or at the local level.

Obtain CARD data:

   .. code-block:: sh
   
      wget https://card.mcmaster.ca/latest/data
      tar -xvf data ./card.json

Local or working directory:

   .. code-block:: sh
   
      rgi load --card_json /path/to/card.json --local

System wide:

   .. code-block:: sh

      rgi load --card_json /path/to/card.json

Check Database Version
-----------------------

Local or working directory:

   .. code-block:: sh
   
      rgi database --version --local

System wide :

   .. code-block:: sh

      rgi database --version
      
Clean Previous or Old Databases
--------------------------------

Local or working directory:

   .. code-block:: sh

      rgi clean --local

System wide:

   .. code-block:: sh 
   
      rgi clean      

RGI main Usage for Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes
------------------------------------------------------------------------------------------------------

.. code-block:: sh

   rgi main -h

.. code-block:: sh

          usage: rgi main [-h] -i INPUT_SEQUENCE -o OUTPUT_FILE [-t {contig,protein}]
                          [-a {DIAMOND,BLAST}] [-n THREADS] [--include_loose] [--local]
                          [--clean] [--debug] [--low_quality]
                          [-d {wgs,plasmid,chromosome,NA}] [-v] [--split_prodigal_jobs]
          
          Resistance Gene Identifier - 4.2.2 - Main
          
          optional arguments:
            -h, --help            show this help message and exit
            -i INPUT_SEQUENCE, --input_sequence INPUT_SEQUENCE
                                  input file must be in either FASTA (contig and
                                  protein) or gzip format! e.g myFile.fasta,
                                  myFasta.fasta.gz
            -o OUTPUT_FILE, --output_file OUTPUT_FILE
                                  output folder and base filename
            -t {contig,protein}, --input_type {contig,protein}
                                  specify data input type (default = contig)
            -a {DIAMOND,BLAST}, --alignment_tool {DIAMOND,BLAST}
                                  specify alignment tool (default = BLAST)
            -n THREADS, --num_threads THREADS
                                  number of threads (CPUs) to use in the BLAST search
                                  (default=32)
            --include_loose       include loose hits in addition to strict and perfect
                                  hits
            --local               use local database (default: uses database in
                                  executable directory)
            --clean               removes temporary files
            --debug               debug mode
            --low_quality         use for short contigs to predict partial genes
            -d {wgs,plasmid,chromosome,NA}, --data {wgs,plasmid,chromosome,NA}
                                  specify a data-type (default = NA)
            -v, --version         prints software version number
            --split_prodigal_jobs
                                  run multiple prodigal jobs simultaneously for contigs
                                  in a fasta file

Running RGI main with Genome or Assembly DNA Sequences
--------------------------------------------------------

Examples use local database, exclude "--local" flag to use a system wide reference database.

Generate Perfect or Strict hits for a genome assembly or genome sequence:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig --local 
      
Include Loose hits:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig --local --include_loose

Short or low quality contigs with partial gene prediction, including Loose hits:

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig --local --low_quality --include_loose

High-performance (e.g. 40 processors) generation of Perfect and Strict hits for high quality genome assembly contigs:

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig --local -a DIAMOND -n 40 --split_prodigal_jobs

Running RGI main with Protein Sequences
--------------------------------------------------------

Examples use local database, exclude "--local" flag to use a system wide reference database.

Generate Perfect or Strict hits for a set of protein sequences:

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/protein_input.fasta --output_file /path/to/output_file --input_type protein --local 

Include Loose hits:

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/protein_input.fasta --output_file /path/to/output_file --input_type protein --local --include_loose

High-performance (e.g. 40 processors) generation of Perfect and Strict hits:

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/protein_input.fasta --output_file /path/to/output_file --input_type protein --local -a DIAMOND -n 40

Running RGI main using GNU Parallel
--------------------------------------------

System wide and writing log files for each input file. Note: add code below to script.sh then run with `./script.sh /path/to/input_files`.

   .. code-block:: sh

      #!/bin/bash
      DIR=`find . -mindepth 1 -type d`
      for D in $DIR; do
            NAME=$(basename $D);
            parallel --no-notice --progress -j+0 'rgi main -i {} -o {.} -n 16 -a diamond --clean --debug > {.}.log 2>&1' ::: $NAME/*.{fa,fasta};
      done

Generating Heat Maps of RGI main Results
------------------------------------------------

.. code-block:: sh

   rgi heatmap -h

.. code-block:: sh

         usage: rgi heatmap [-h] -i INPUT
                            [-cat {drug_class,resistance_mechanism,gene_family}] [-f]
                            [-o OUTPUT] [-clus {samples,genes,both}]
                            [-d {plain,fill,text}] [--debug]
         
         Creates a heatmap when given multiple RGI results.
         
         optional arguments:
           -h, --help            show this help message and exit
           -i INPUT, --input INPUT
                                 Directory containing the RGI .json files (REQUIRED)
           -cat {drug_class,resistance_mechanism,gene_family}, --category {drug_class,resistance_mechanism,gene_family}
                                 The option to organize resistance genes based on a
                                 category.
           -f, --frequency       Represent samples based on resistance profile.
           -o OUTPUT, --output OUTPUT
                                 Name for the output EPS and PNG files. The number of
                                 files run will automatically be appended to the end of
                                 the file name. (default=RGI_heatmap)
           -clus {samples,genes,both}, --cluster {samples,genes,both}
                                 Option to use SciPy's hiearchical clustering algorithm
                                 to cluster rows (AMR genes) or columns (samples).
           -d {plain,fill,text}, --display {plain,fill,text}
                                 Specify display options for categories
                                 (deafult=plain).
           --debug               debug mode


RGI heatmap produces EPS and PNG image files.

Generate a heat map from pre-compiled RGI main JSON files, samples and AMR genes organized alphabetically:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --output /path/to/output_file
            
Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by AMR gene family:            

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --output /path/to/output_file -cat gene_family -clus samples

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by Drug Class:            

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --output /path/to/output_file -cat drug_class -clus samples

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by distribution among samples:            

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --output /path/to/output_file -clus both
            
Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome (with histogram used for abundance of identical resistomes) and AMR genes organized by distribution among samples:            

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --output /path/to/output_file -clus both -f

Run RGI from Docker
-------------------

First you you must either pull the Docker container from dockerhub (latest CARD version automatically installed):

  .. code-block:: sh

        docker pull finlaymaguire/rgi

Or alternatively, build it locally from the Dockerfile (latest CARD version automatically installed):

  .. code-block:: sh

        git clone https://github.com/arpcard/rgi
        docker build -t arpcard/rgi rgi

Then you can either run interactively (mounting a local directory called `rgi_data` in your current directory to `/data/` within the container:

  .. code-block:: sh

        docker run -i -v $PWD/rgi_data:/data -t arpcard/rgi bash

Or you can directly run the container as an executable with `$RGI_ARGS` being any of the commands described above. Remember paths to input and outputs files are relative to the container (i.e. `/data/` if mounted as above):

  .. code-block:: sh
        
        docker run -v $PWD/rgi_data:/data arpcard/rgi $RGI_ARGS

Install RGI from Conda
-------------------

Search for RGI package and show available versions:

  .. code-block:: sh
        
        $ conda search --channel bioconda rgi

Install RGI package:

  .. code-block:: sh
        
        $ conda install --channel bioconda rgi

Install RGI specific version:

  .. code-block:: sh
        
        $ conda install --channel bioconda rgi=3.1.1

Remove RGI package:

  .. code-block:: sh
        
        $ conda remove --channel bioconda rgi

Overview of Tab-Delimited Output
-----------------------------------


+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    ORF_ID                                                | Open Reading Frame identifier (internal to RGI)|
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Contig                                                | Source Sequence                                |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Start                                                 | Start co-ordinate of ORF                       |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Stop                                                  | End co-ordinate of ORF                         |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Orientation                                           | Strand of ORF                                  |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Cut_Off                                               | RGI Detection Paradigm                         |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Pass_Bitscore                                         | STRICT detection model bitscore value cut-off  |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Best_Hit_Bitscore                                     | Bitscore value of match to top hit in CARD     |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Best_Hit_ARO                                          | ARO term of top hit in CARD                    |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Best_Identities                                       | Percent identity of match to top hit in CARD   |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    ARO                                                   | ARO accession of top hit in CARD               |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Model_type                                            | CARD detection model type                      |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|                                                          | Mutations observed in the ARO term of top hit  |
|    SNPs_in_Best_Hit_ARO                                  | in CARD (if applicable)                        |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|                                                          | Mutations observed in ARO terms of other hits  |
|    Other_SNPs                                            | indicated by model id (if applicable)          |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Drug Class                                            | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Resistance Mechanism                                  | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    AMR Gene Family                                       | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Predicted_DNA                                         | ORF predicted nucleotide sequence              |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Predicted_Protein                                     | ORF predicted protein sequence                 |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    CARD_Protein_Sequence                                 | Protein sequence of top hit in CARD            |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       | Calculated as percentage                       |
|                                                          | (length of ORF protein /                       |
|    Percentage Length of Reference Sequence               | length of CARD reference protein)              |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    ID                                                    | HSP identifier (internal to RGI)               |
+----------------------------------------------------------+------------------------------------------------+
| ::                                                       |                                                |
|    Model_id                                              | CARD detection model id                        |
+----------------------------------------------------------+------------------------------------------------+

