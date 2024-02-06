K-mer Prediction of Pathogen-of-Origin for AMR Genes
----------------------------------------------------

CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ and `Prevalence Data <https://card.mcmaster.ca/prevalence>`_ (see above) provides a data set of AMR alleles and their distribution among pathogens and plasmids. CARD's k-mer classifiers sub-sample these sequences to identify k-mers (default length 61 bp) that are uniquely found within AMR alleles of individual pathogen species, pathogen genera, pathogen-restricted plasmids, or promiscuous plasmids. CARD's k-mer classifiers can then be used to predict pathogen-of-origin for matches found by RGI for genomes, genome assemblies, metagenomic contigs, or metagenomic reads.

**CARD's k-mer classifiers assume the data submitted for analysis has been predicted to encode AMR genes, via RGI or another AMR bioinformatic tool. The k-mer data set was generated from and is intended exclusively for AMR sequence space.** As above, the reported results are entirely dependant upon the curated AMR detection models in CARD, the algorithms available in RGI, and the pathogens & sequences sampled during generation of CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ and `Prevalence Data <https://card.mcmaster.ca/prevalence>`_.

Using RGI kmer_query
--------------------

**This is an unpublished algorithm undergoing beta-testing.**

As outlined above, CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ and `Prevalence Data <https://card.mcmaster.ca/prevalence>`_ provide a data set of AMR alleles and their distribution among pathogens and plasmids. CARD's k-mer classifiers sub-sample these sequences to identify k-mers that are uniquely found within AMR alleles of individual pathogen species, pathogen genera, pathogen-restricted plasmids, or promiscuous plasmids. The default k-mer length is 61 bp (based on unpublished analyses), available as downloadable, pre-compiled k-mer sets at the CARD website.

CARD's k-mer classifiers assume the data submitted for analysis has been predicted to encode AMR genes, via RGI or another AMR bioinformatic tool. The k-mer data set was generated from and is intended exclusively for AMR sequence space. To be considered for a taxonomic prediction, individual sequences (e.g. FASTA, RGI predicted ORF, metagenomic read) must pass the *--minimum* coverage value (default of 10, i.e. the number of k-mers in a sequence that need to match a single category, for both taxonomic and genomic classifications, in order for a classification to be made for that sequence). Subsequent classification is based on the following logic tree:

.. image:: /images/kmerlogic.jpg

.. code-block:: sh

   rgi kmer_query -h

.. code-block:: sh

				usage: rgi kmer_query [-h] -i INPUT [--bwt] [--rgi] [--fasta] -k K [-m MIN]
				                      [-n THREADS] -o OUTPUT [--local] [--debug]

				Resistance Gene Identifier - 6.0.2 - Kmer Query

				Tests sequenes using CARD*kmers

				optional arguments:
				  -h, --help            show this help message and exit
				  -i INPUT, --input INPUT
				                        Input file (bam file from RGI*BWT, json file of RGI results, fasta file of sequences)
				  --bwt                 Specify if the input file for analysis is a bam file generated from RGI*BWT
				  --rgi                 Specify if the input file is a RGI results json file
				  --fasta               Specify if the input file is a fasta file of sequences
				  -k K, --kmer_size K   length of k
				  -m MIN, --minimum MIN
				                        Minimum number of kmers in the called category for the classification to be made (default=10).
				  -n THREADS, --threads THREADS
				                        number of threads (CPUs) to use (default=1)
				  -o OUTPUT, --output OUTPUT
				                        Output file name.
				  --local               use local database (default: uses database in executable directory)
				  --debug               debug mode

If you have not already done so, you must load CARD reference data for these commands to work. First, remove any previous loads:

   .. code-block:: sh

      rgi clean --local

Download CARD data:

   .. code-block:: sh

      wget https://card.mcmaster.ca/latest/data
      tar -xvf data ./card.json

Load into local or working directory:

   .. code-block:: sh

      rgi load --card_json /path/to/card.json --local

Also pre-process these reference data for metagenomics reads (note that the filename *card_database_v3.0.1.fasta* depends on the version of CARD data downloaded, please adjust accordingly):

   .. code-block:: sh

      rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1
      rgi load -i /path/to/card.json --card_annotation card_database_v3.0.1.fasta --local

The pre-compiled 61 bp k-mers are available via CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_:

   .. code-block:: sh

      wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
      mkdir -p wildcard
      tar -xjf wildcard_data.tar.bz2 -C wildcard
      gunzip wildcard/*.gz

Load k-mers:

   .. code-block:: sh

      rgi load --card_json /path/to/card.json
        --kmer_database /path/to/wildcard/61_kmer_db.json
        --amr_kmers /path/to/wildcard/all_amr_61mers.txt --kmer_size 61
        --local --debug > kmer_load.61.log 2>&1

CARD k-mer Classifier analysis of an individual FASTA file (e.g. using 8 processors, minimum k-mer coverage of 10):

.. code-block:: sh

   rgi kmer_query --fasta --kmer_size 61 --threads 8 --minimum 10
    --input /path/to/nucleotide_input.fasta --output /path/to/output_file --local

CARD k-mer Classifier analysis of Genome or Assembly DNA Sequences RGI main results (e.g. using 8 processors, minimum k-mer coverage of 10):

.. code-block:: sh

   rgi kmer_query --rgi --kmer_size 61 --threads 8 --minimum 10
    --input /path/to/rgi_main.json --output /path/to/output_file --local

CARD k-mer Classifier analysis of Metagenomics RGI btw results (e.g. using 8 processors, minimum k-mer coverage of 10):

.. code-block:: sh

   rgi kmer_query --bwt --kmer_size 61 --threads 8 --minimum 10
    --input /path/to/rgi_bwt.bam --output /path/to/output_file --local

CARD k-mer Classifier Output
````````````````````````````

CARD k-mer classifier output differs between genome/gene and metagenomic data:

CARD k-mer Classifier Output for a FASTA file
`````````````````````````````````````````````

+----------------------------------------------------------+----------------------------------------------------+
|    Field                                                 | Contents                                           |
+==========================================================+====================================================+
|    Sequence                                              | Sequence defline in the FASTA file                 |
+----------------------------------------------------------+----------------------------------------------------+
|    Total # kmers                                         | Total # k-mers in the sequence                     |
+----------------------------------------------------------+----------------------------------------------------+
|    # of AMR kmers                                        | Total # AMR k-mers in the sequence                 |
+----------------------------------------------------------+----------------------------------------------------+
|    CARD kmer Prediction                                  | Taxonomic prediction, with indication if the k-mers|
|                                                          | are known exclusively from chromosomes, exclusively|
|                                                          | from plasmids, or can be found in either           |
|                                                          | chromosomes or plasmids                            |
+----------------------------------------------------------+----------------------------------------------------+
|    Taxonomic kmers                                       | Number of k-mer hits broken down by taxonomy       |
+----------------------------------------------------------+----------------------------------------------------+
|    Genomic kmers                                         | Number of k-mer hits exclusive to chromosomes,     |
|                                                          | exclusively to plasmids, or found in either        |
|                                                          | chromosomes or plasmids                            |
+----------------------------------------------------------+----------------------------------------------------+

CARD k-mer Classifier Output for RGI main results
`````````````````````````````````````````````````

+----------------------------------------------------------+----------------------------------------------------+
|    Field                                                 | Contents                                           |
+==========================================================+====================================================+
|    ORF_ID                                                | Open Reading Frame identifier (from RGI results)   |
+----------------------------------------------------------+----------------------------------------------------+
|    Contig                                                | Source Sequence (from RGI results)                 |
+----------------------------------------------------------+----------------------------------------------------+
|    Cut_Off                                               | RGI Detection Paradigm (from RGI results)          |
+----------------------------------------------------------+----------------------------------------------------+
|    CARD kmer Prediction                                  | Taxonomic prediction, with indication if the k-mers|
|                                                          | are known exclusively from chromosomes, exclusively|
|                                                          | from plasmids, or can be found in either           |
|                                                          | chromosomes or plasmids                            |
+----------------------------------------------------------+----------------------------------------------------+
|    Taxonomic kmers                                       | Number of k-mer hits broken down by taxonomy       |
+----------------------------------------------------------+----------------------------------------------------+
|    Genomic kmers                                         | Number of k-mer hits exclusive to chromosomes,     |
|                                                          | exclusively to plasmids, or found in either        |
|                                                          | chromosomes or plasmids                            |
+----------------------------------------------------------+----------------------------------------------------+

CARD k-mer Classifier Output for RGI bwt results
````````````````````````````````````````````````

As with RGI bwt analysis, output is produced at both the allele and gene level:

+----------------------------------------------------------+----------------------------------------------------+
|    Field                                                 | Contents                                           |
+==========================================================+====================================================+
|    Reference Sequence / ARO term                         | Reference allele or gene ARO term to which reads   |
|                                                          | have been mapped                                   |
+----------------------------------------------------------+----------------------------------------------------+
|    Mapped reads with kmer DB hits                        | **Number of reads** classified                     |
+----------------------------------------------------------+----------------------------------------------------+
|    CARD kmer Prediction                                  | **Number of reads** classified for each allele or  |
|                                                          | gene, with indication if the k-mers are known      |
|                                                          | exclusively from chromosomes, exclusively from     |
|                                                          | plasmids, or can be found in either                |
+----------------------------------------------------------+----------------------------------------------------+
|    Subsequent fields                                     | Detected k-mers within the context of the k-mer    |
|                                                          | logic tree                                         |
+----------------------------------------------------------+----------------------------------------------------+

Building Custom k-mer Classifiers
`````````````````````````````````

**This is an unpublished algorithm undergoing beta-testing.**

You must `Load CARD Reference Data`_ for these commands to work.

As outlined above, CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ and `Prevalence Data <https://card.mcmaster.ca/prevalence>`_ provide a data set of AMR alleles and their distribution among pathogens and plasmids. CARD's k-mer classifiers sub-sample these sequences to identify k-mers that are uniquely found within AMR alleles of individual pathogen species, pathogen genera, pathogen-restricted plasmids, or promiscuous plasmids. The default k-mer length is 61 bp (based on unpublished analyses), available as downloadable, pre-compiled k-mer sets at the CARD website, but users can also use RGI to create k-mers of any length. **Warning**: this is computationally intensive.

.. code-block:: sh

   rgi kmer_build -h

.. code-block:: sh

				usage: rgi kmer_build [-h] [-i INPUT_DIRECTORY] -c CARD_FASTA -k K [--skip]
				                      [-n THREADS] [--batch_size BATCH_SIZE]

				Resistance Gene Identifier - 6.0.2 - Kmer Build

				Builds the kmer sets for CARD*kmers

				optional arguments:
				  -h, --help            show this help message and exit
				  -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
				                        input directory of prevalence data
				  -c CARD_FASTA, --card CARD_FASTA
				                        fasta file of CARD reference sequences. If missing, run 'rgi card_annotation' to generate.
				  -k K                  k-mer size (e.g., 61)
				  --skip                skips the concatenation and splitting of the CARD*R*V sequences.
				  -n THREADS, --threads THREADS
				                        number of threads (CPUs) to use (default=1)
				  --batch_size BATCH_SIZE
				                        number of kmers to query at a time using pyahocorasick--the greater the number the more memory usage (default=100,000)

Example generation of 31 bp k-mers using 20 processors (note that the filename *card_database_v3.0.1.fasta* depends on the version of CARD data downloaded, please adjust accordingly):

.. code-block:: sh

   rgi kmer_build --input_directory /path/to/wildcard
    --card card_database_v3.0.1.fasta -k 31 --threads 20 --batch_size 100000

The *--skip* flag can be used if you are making k-mers a second time (33 bp in the example below) to avoid re-generating intermediate files (note that the filename *card_database_v3.0.1.fasta* depends on the version of CARD data downloaded, please adjust accordingly):

.. code-block:: sh

   rgi kmer_build --input_directory /path/to/wildcard
    --card card_database_v3.0.1.fasta -k 33 --threads 20 --batch_size 100000 --skip
    