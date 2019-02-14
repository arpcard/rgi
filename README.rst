Resistance Gene Identifier (RGI) 
--------------------------------------------

This application is used to predict resistome(s) from protein or nucleotide data based on homology and SNP models. The application uses data from `CARD database <https://card.mcmaster.ca/>`_.

.. image:: https://travis-ci.org/arpcard/rgi.svg?branch=master
    :target: https://travis-ci.org/arpcard/rgi

Table of Contents
-------------------------------------

- `License`_
- `Requirements`_
- `Install dependencies`_
- `Install RGI from project root`_
- `Running RGI tests`_
- `Help menu`_
- `Usage`_
- `Load card.json`_
- `Check database version`_
- `Run RGI`_
- `Run RGI using GNU parallel`_
- `Running RGI with short contigs to predict partial genes`_
- `Clean previous or old databases`_
- `RGI Heatmap`_
- `Run RGI from docker`_
- `Tab-delimited results file`_
- `Support & Bug Reports`_

License
--------
Use or reproduction of these materials, in whole or in part, by any non-academic organization whether or not for non-commercial (including research) or commercial purposes is prohibited, except with written permission of McMaster University. Commercial uses are offered only pursuant to a written license and user fee. To obtain permission and begin the licensing process, see `CARD website <https://card.mcmaster.ca/about>`_

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

Install dependencies
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

Install RGI from project root
-----------------------------

.. code-block:: sh

   pip3 install .

or

.. code-block:: sh

   python3 setup.py build
   python3 setup.py test
   python3 setup.py install

Running RGI tests
-------------------
.. code-block:: sh
   
   cd tests
   pytest -v -rxs

Help menu
-------------------

.. code-block:: sh

   rgi --help

Usage
-------------------

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


Load card.json 
-------------------

- local or working directory

   .. code-block:: sh
   
      rgi load --card_json /path/to/card.json --local

- system wide 

   .. code-block:: sh

      rgi load --card_json /path/to/card.json

Check database version
-----------------------

- local or working directory

   .. code-block:: sh
   
      rgi database --version --local

- system wide 

   .. code-block:: sh

      rgi database --version

Run RGI 
----------------------

- local or working directory

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/protein_input.fasta --output_file /path/to/output_file --input_type protein --local 

- system wide 

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig

Run RGI using GNU parallel
--------------------------------------------

- system wide and writing log files for each input file. (Note add code below to script.sh then run with `./script.sh /path/to/input_files`)

   .. code-block:: sh

      #!/bin/bash
      DIR=`find . -mindepth 1 -type d`
      for D in $DIR; do
            NAME=$(basename $D);
            parallel --no-notice --progress -j+0 'rgi main -i {} -o {.} -n 16 -a diamond --clean --debug > {.}.log 2>&1' ::: $NAME/*.{fa,fasta};
      done



Running RGI with short contigs to predict partial genes 
--------------------------------------------------------

- local or working directory

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --local --low_quality 

- system wide 

   .. code-block:: sh
   
      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --low_quality


Clean previous or old databases
--------------------------------

- local or working directory

   .. code-block:: sh

      rgi clean --local

- system wide 

   .. code-block:: sh 
   
      rgi clean      

RGI Heatmap
------------

- Default Heatmap

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
       
- Heatmap with `AMR Gene Family` categorization

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --category gene_family

- Heatmap with `AMR Gene Family` categorization and fill display

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --category gene_family --display fill

- Heatmap with `AMR Gene Family` categorization and coloured y-axis labels display

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --category gene_family --display text


- Heatmap with frequency display enabled

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --frequency

- Heatmap with drug class category and frequency enabled

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --category drug_class --frequency --display text

- Heatmap with samples and genes clustered

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --cluster both

- Heatmap with resistance mechanism categorization and clustered samples

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/ --cluster samples --category resistance_mechanism --display fill


Run RGI from docker
-------------------

- First you you must either pull the docker container from dockerhub (latest CARD version automatically installed)

  .. code-block:: sh

        docker pull finlaymaguire/rgi

- Or Alternatively, build it locally from the Dockerfile (latest CARD version automatically installed)

  .. code-block:: sh

        git clone https://github.com/arpcard/rgi
        docker build -t arpcard/rgi rgi

- Then you can either run interactively (mounting a local directory called `rgi_data` in your current directory
  to `/data/` within the container

  .. code-block:: sh

        docker run -i -v $PWD/rgi_data:/data -t arpcard/rgi bash

- Or you can directly run the container as an executable with `$RGI_ARGS` being any of the commands described above. Remember paths to input and outputs files are relative to the container (i.e. `/data/` if mounted as above).

  .. code-block:: sh
        
        docker run -v $PWD/rgi_data:/data arpcard/rgi $RGI_ARGS
       
Tab-delimited results file
---------------------------

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



Support & Bug Reports
----------------------

Please log an issue on `github issue <https://github.com/arpcard/rgi/issues>`_.

You can email the CARD curators or developers directly at `card@mcmaster.ca <mailto:card@mcmaster.ca>`_, via Twitter at `@arpcard <http://www.twitter.com/arpcard>`_.

