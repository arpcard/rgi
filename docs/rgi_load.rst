.. contents:: Page Contents
   :local:

Loading CARD Reference Data
````````````````````````````

.. code-block:: sh

				usage: rgi load [-h] -i CARD_JSON [--card_annotation CARD_ANNOTATION]
				                [--card_annotation_all_models CARD_ANNOTATION_ALL_MODELS]
				                [--wildcard_annotation WILDCARD_ANNOTATION]
				                [--wildcard_annotation_all_models WILDCARD_ANNOTATION_ALL_MODELS]
				                [--wildcard_index WILDCARD_INDEX]
				                [--wildcard_version WILDCARD_VERSION]
				                [--baits_annotation BAITS_ANNOTATION]
				                [--baits_index BAITS_INDEX] [--kmer_database KMER_DATABASE]
				                [--amr_kmers AMR_KMERS] [--kmer_size KMER_SIZE] [--local]
				                [--debug] [--include_other_models]

				Resistance Gene Identifier - 6.0.2 - Load

				optional arguments:
				  -h, --help            show this help message and exit
				  -i CARD_JSON, --card_json CARD_JSON
				                        must be a card database json file
				  --card_annotation CARD_ANNOTATION
				                        annotated reference FASTA for protein homolog models
				                        only, created using rgi card_annotation
				  --card_annotation_all_models CARD_ANNOTATION_ALL_MODELS
				                        annotated reference FASTA which includes all models
				                        created using rgi card_annotation
				  --wildcard_annotation WILDCARD_ANNOTATION
				                        annotated reference FASTA for protein homolog models
				                        only, created using rgi wildcard_annotation
				  --wildcard_annotation_all_models WILDCARD_ANNOTATION_ALL_MODELS
				                        annotated reference FASTA which includes all models
				                        created using rgi wildcard_annotation
				  --wildcard_index WILDCARD_INDEX
				                        wildcard index file (index-for-model-sequences.txt)
				  --wildcard_version WILDCARD_VERSION
				                        specify variants version used
				  --baits_annotation BAITS_ANNOTATION
				                        annotated reference FASTA
				  --baits_index BAITS_INDEX
				                        baits index file (baits-probes-with-sequence-info.txt)
				  --kmer_database KMER_DATABASE
				                        json of kmer database
				  --amr_kmers AMR_KMERS
				                        txt file of all amr kmers
				  --kmer_size KMER_SIZE
				                        kmer size if loading kmer files
				  --local               use local database (default: uses database in
				                        executable directory)
				  --debug               debug mode

Depending upon the type of analysis you wish to perform, different sets of CARD reference data first need to be loaded into RGI. By default, these data will be loaded at the system-wide level, i.e. available to all users alongside a system-wide RGI installation, but they can alternatively be loaded for the local user directory using the --local flag. Steps for loading required data are outlined below in sections describing different types of analysis (all using --local in their examples), but below are examples of loading the canonical CARD reference data either system-wide or locally.

First download the latest AMR reference data from CARD:

   .. code-block:: sh

      wget https://card.mcmaster.ca/latest/data
      tar -xvf data ./card.json

Load in Local or working directory:

   .. code-block:: sh

      rgi load --card_json /path/to/card.json --local

Load System wide:

   .. code-block:: sh

      rgi load --card_json /path/to/card.json

Check Database Version
``````````````````````

Local or working directory:

   .. code-block:: sh

      rgi database --version --local

System wide :

   .. code-block:: sh

      rgi database --version

Clean Previous or Old Databases
````````````````````````````````

Local or working directory:

   .. code-block:: sh

      rgi clean --local

System wide:

   .. code-block:: sh

      rgi clean

Bulk Load All Reference Data
`````````````````````````````

The examples in this documentation outline best practices for loading of CARD reference data for each possible type of analysis. If you wish to bulk load all possible CARD reference data to allow on-the-fly switching between different types of analysis, here are all of the steps combined:

Remove any previous loads:

   .. code-block:: sh

      rgi clean --local

Download CARD and WildCARD data:

   .. code-block:: sh

      wget https://card.mcmaster.ca/latest/data
      tar -xvf data ./card.json
      wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
      mkdir -p wildcard
      tar -xjf wildcard_data.tar.bz2 -C wildcard
      gunzip wildcard/*.gz

Create annotation files (note that the parameter *version_number* depends upon the versions of WildCARD data downloaded, please adjust accordingly):

   .. code-block:: sh

      rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1
      rgi wildcard_annotation -i wildcard --card_json /path/to/card.json
        -v version_number > wildcard_annotation.log 2>&1

Load all data into RGI (note that the FASTA filenames plus the parameter *version_number* depend on the versions of CARD and WildCARD data downloaded, please adjust accordingly):

   .. code-block:: sh

     rgi load \
       --card_json /path/to/card.json \
       --debug --local \
       --card_annotation card_database_v3.2.4.fasta \
       --card_annotation_all_models card_database_v3.2.4_all.fasta \
       --wildcard_annotation wildcard_database_v4.0.0.fasta \
       --wildcard_annotation_all_models wildcard_database_v4.0.0_all.fasta \
       --wildcard_index /path/to/wildcard/index-for-model-sequences.txt \
       --wildcard_version 4.0.0 \
       --amr_kmers /path/to/wildcard/all_amr_61mers.txt \
       --kmer_database /path/to/wildcard/61_kmer_db.json \
       --kmer_size 61


