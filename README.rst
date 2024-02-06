|build-status| |docs|

.. |build-status| image:: https://travis-ci.org/arpcard/rgi.svg?branch=master
    :alt: build status
    :scale: 100%
    :target: https://travis-ci.org/arpcard/rgi

.. |build-status| image:: https://github.com/arpcard/rgi/actions/workflows/build.yml/badge.svg?branch=master
		 :alt: Workflow status badge
		 :scale: 100%
		 :target: https://github.com/arpcard/rgi/actions/workflows/build.yml

.. |docs| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :alt: Documentation
    :scale: 100%
    :target: http://bioconda.github.io/recipes/rgi/README.html

====================================
The Resistance Gene Identifier (RGI)
====================================

This application is used to predict antibiotic resistome(s) from protein or nucleotide data based on homology and SNP models. The application uses reference data from the `Comprehensive Antibiotic Resistance Database (CARD) <https://card.mcmaster.ca/>`_.

RGI analyses can be performed via the CARD website `RGI portal <https://card.mcmaster.ca/analyze/rgi>`_, via use of a `Galaxy wrapper <https://toolshed.g2.bx.psu.edu/view/card/rgi/715bc9aeef69>`_ for the `Galaxy <https://galaxyproject.org/tutorials/g101>`_ platform, or alternatively you can install RGI from Conda or run RGI from Docker (see below). The instructions below discuss use of RGI at the command line, following a general overview of how RGI works for genomes, genome assemblies, proteomes, and metagenomic sequencing.

May 2023: Chan Zuckerberg ID (CZ ID) has implemented a web-based platform for RGI analysis of assembled contigs (FASTA) or metagenomic sequencing reads (FASTQ): `CZ ID AMR Pipeline Workflow <https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow>`_.

**CARD reference sequences and significance cut-offs are under constant curation - as CARD curation evolves, the results of RGI evolve.**

* `CARD Frequency Asked Questions <https://github.com/arpcard/FAQ>`_
* YouTube Video Lecture - Canadian Bioinformatics Workshops 2023: `Infectious Disease Genomic Epidemiology - Antimicrobial Resistant Gene (AMR) Analysis using CARD & RGI <https://www.youtube.com/watch?v=FvOCDlcYaTo&list=PL3izGL6oi0S8RG8vnwLXFznzJnKh8OR8F&index=6>`_

.. contents::

Overview of RGI
===============

* `Analyzing Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes </docs/rgi_main.rst>`_ (a.k.a. RGI main)
* `Analyzing Metagenomic Reads </docs/rgi_bwt.rst>`_ (a.k.a. RGI bwt)
* `K-mer Prediction of Pathogen-of-Origin for AMR Genes </docs/rgi_kmer.rst>`_ (beta-testing)

License
--------

Use or reproduction of these materials, in whole or in part, by any commercial organization whether or not for non-commercial (including research) or commercial purposes is prohibited, except with written permission of McMaster University. Commercial uses are offered only pursuant to a written license and user fee. To obtain permission and begin the licensing process, see the `CARD website <https://card.mcmaster.ca/about>`_.

Citation
--------

Alcock et al. 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690-D699 [`PMID 36263822 <https://www.ncbi.nlm.nih.gov/pubmed/36263822>`_]

Support & Bug Reports
----------------------

Please log an issue on `github issue <https://github.com/arpcard/rgi/issues>`_.

You can email the CARD curators or developers directly at `card@mcmaster.ca <mailto:card@mcmaster.ca>`_.

---------------------


Installation
============

Recommended installation method for most users is via Conda or Docker.
This will handle dependency management and ensure installation of the
correct version of RGI's external dependencies e.g., BLAST, DIAMOND.

Install RGI from Conda
----------------------

Install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_ on your system if not already available.

Install `mamba` from `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ on your system if not already available.

Search for RGI package and show available versions:

  .. code-block:: sh

        $ mamba search --channel conda-forge --channel bioconda --channel defaults rgi

Create a new Conda environment

  .. code-block:: sh

        $ mamba create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi

Install RGI package:

  .. code-block:: sh

        $ mamba install --channel conda-forge --channel bioconda --channel defaults rgi

Install RGI specific version:

  .. code-block:: sh

        $ mamba install --channel conda-forge --channel bioconda --channel defaults rgi=5.1.1

Remove RGI package:

  .. code-block:: sh

        $ mamba remove rgi


Install RGI using Docker/Singularity
------------------------------------

RGI is available via biocontainers full installed with all
databases appropriately loaded.

Install `docker <https://docs.docker.com/get-docker/>`_ on your system if not already available

- Pull the Docker container from biocontainers (built from Conda package at https://quay.io/repository/biocontainers/rgi?tab=tags&tag=latest).

    .. code-block:: sh

        docker pull quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0

- RGI can be executed from the container as follows:

    .. code-block:: sh

        docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0 rgi -h


Install Development Version
---------------------------

Install Dependencies
````````````
The following conda command will install all RGI dependencies (listed below):

.. code-block:: sh

    git clone https://github.com/arpcard/rgi
    conda env create -f conda_env.yml
    conda activate rgi


- `Python 3.6 <https://www.python.org/>`_
- `NCBI BLAST 2.14.0 <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
- `zlib <https://bitbucket.org/gutworth/six>`_
- `Prodigal 2.6.3 <https://github.com/hyattpd/prodigal/wiki/Installation>`_
- `DIAMOND 0.8.36 <https://github.com/bbuchfink/diamond>`_
- `Biopython 1.78 <https://biopython.org/>`_
- `filetype 1.0.0+ <https://pypi.org/project/filetype/>`_
- `pytest 3.0.0+ <https://docs.pytest.org/en/latest/>`_
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
- `KMA 1.3.4 <https://bitbucket.org/genomicepidemiology/kma/src/master>`_


Install RGI
```````````

.. code-block:: sh

   pip install git+https://github.com/arpcard/rgi.git

or

.. code-block:: sh

   python setup.py build
   python setup.py test
   python setup.py install

Running RGI Tests
`````````````````
.. code-block:: sh

   cd tests
   pytest -v -rxs

-------------------

RGI Usage Documentation
=======================

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
               auto_load Automatically loads CARD database, annotations and k-mer database
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
               bwt                   Align reads to CARD and in silico predicted allelic variants (beta)

               ---------------------------------------------------------------------------------------
               Baits validation
               ---------------------------------------------------------------------------------------
               tm                    Baits Melting Temperature

               ---------------------------------------------------------------------------------------
               Annotations
               ---------------------------------------------------------------------------------------
               card_annotation       Create fasta files with annotations from card.json
               wildcard_annotation   Create fasta files with annotations from variants
               baits_annotation      Create fasta files with annotations from baits (experimental)
               remove_duplicates     Removes duplicate sequences (experimental)

               ---------------------------------------------------------------------------------------
               Pathogen of origin
               ---------------------------------------------------------------------------------------

               kmer_build            Build AMR specific k-mers database used for pathogen of origin (beta)
               kmer_query            Query sequences against AMR k-mers database to predict pathogen of origin (beta)

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


RGI Databases
--------------

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

Running RGI on Compute Canada Serial Farm
`````````````````````````````````````````

**Order of operations**

.. code-block:: sh

   ## Running jobs on computecanada using serial farm method

   - `rgi bwt` was used as example.

   ### step 1:

   - update make_table_dat.sh to construct arguments for commands

   ### step 2:

   - update eval command in job_script.sh to match your tool and also load appropriate modules

   ### step 3:

   - create table.dat using script make_table_dat.sh with inputs files in all_samples directory
   ./make_table_dat.sh ./all_samples/ > table.dat

   ### step 4:

   - submit multiple jobs using for_loop.sh

   ### Resource:

   - https://docs.computecanada.ca/wiki/Running_jobs#Serial_job


**Update the make_table_dat.sh**

.. code-block:: sh

   DIR=`find . -mindepth 1 -type d`
   for D in $DIR; do
         directory=$(basename $D);
         for file in $directory/*; do
           filename=$(basename $file);
         if [[ $filename = *"_pass_1.fastq.gz"* ]]; then
               read1=$(basename $filename);
                base=(${read1//_pass_1.fastq.gz/ });
                #echo "--read_one $(pwd)/$directory/${base}_pass_1.fastq.gz --read_two $(pwd)/$directory/${base}_pass_2.fastq.gz -o $(pwd)/$directory/${base} -n 16 --aligner bowtie2 --debug"
            echo "--read_one $(pwd)/$directory/${base}_pass_1.fastq.gz --read_two $(pwd)/$directory/${base}_pass_2.fastq.gz -o $(pwd)/$directory/${base}_wild -n 8 --aligner bowtie2 --debug --include_wildcard"
         fi
         done
    done

This block of code is used to generate the arguments for serial farming. In this example, rgi bwt is used, however depending on the job you are running you may update it according to your specifications.

**Update the job_script.sh to match used tool**

.. code-block:: sh

   #SBATCH --account=def-mcarthur
   #SBATCH --time=120
   #SBATCH --job-name=rgi_bwt
   #SBATCH --cpus-per-task=8
   #SBATCH --mem-per-cpu=2048M
   #SBATCH --mail-user=raphenar@mcmaster.ca
   #SBATCH --mail-type=ALL

   # Extracing the $I_FOR-th line from file $TABLE:
   LINE=`sed -n ${I_FOR}p "$TABLE"`

   # Echoing the command (optional), with the case number prepended:
   #echo "$I_FOR; $LINE"

   # load modules
   module load nixpkgs/16.09 python/3.6.3 gcc/5.4.0 blast+/2.6.0 prodigal diamond/0.8.36 bowtie2  samtools bamtools bedtools bwa

   # execute command
   #eval "$LINE"
   #echo "rgi bwt $LINE"
   eval "rgi bwt $LINE"

Update this block of code according to which tool you want to use. In this example, rgi bwt is shown, however for your use-case, you may update it accordingly.

**Creating the table.dat**

To create the table.dat, use the script made before named make_table_dat.sh along with the path to the directory containing all your inputs as an argument. Output to table.dat.

.. code-block:: sh

   ./make_table_dat.sh ./all_samples/ > table.dat

**Submit multiple jobs using for_loop.sh**

This script is used once all the previous steps are completed. This script allows you to submit multiple jobs into Compute Canada for rgi.

.. code-block:: sh

   # Simplest case - using for loop to submit a serial farm
   # The input file table.dat contains individual cases - one case per line
   export TABLE=table.dat

   # Total number of cases (= number of jobs to submit):
   N_cases=$(cat "$TABLE" | wc -l)

   # Submitting one job per case using the for loop:
   for ((i=1; i<=$N_cases; i++))
    do
    # Using environment variable I_FOR to communicate the case number to individual jobs:
    export I_FOR=$i
    sbatch job_script.sh
   done

**Resources**

More information on serial farming on Compute Canada can be found here_.

.. _here: https://docs.computecanada.ca/wiki/Running_jobs#Serial_job

