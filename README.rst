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
* `CBW 2024 Infectious Disease Genomic Epidemiology - Antimicrobial Resistant Gene (AMR) Analysis using CARD & RGI <https://www.youtube.com/watch?v=Z4gAAYRExSg&list=PL3izGL6oi0S_e5T8qx-74WRaMR5K5U8V5&index=8>`_

Overview and Use of RGI
=======================

* `Help Menu and Usage </docs/rgi_help.rst>`_
* `Loading CARD Reference Databases </docs/rgi_load.rst>`_
* `Analyzing Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes </docs/rgi_main.rst>`_ (a.k.a. RGI main)
* `Analyzing Metagenomic Reads </docs/rgi_bwt.rst>`_ (a.k.a. RGI bwt)
* `K-mer Prediction of Pathogen-of-Origin for AMR Genes </docs/rgi_kmer.rst>`_ (beta-testing)
* `Running RGI on Digital Alliance (formerly Compute Canada) Serial Farm </docs/rgi_digital_alliance.rst>`_

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

        mamba search --channel conda-forge --channel bioconda --channel defaults rgi

Create a new Conda environment

  .. code-block:: sh

        mamba create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi

Install RGI package:

  .. code-block:: sh

        mamba install --channel conda-forge --channel bioconda --channel defaults rgi

Install RGI specific version:

  .. code-block:: sh

        mamba install --channel conda-forge --channel bioconda --channel defaults rgi=5.1.1

Remove RGI package:

  .. code-block:: sh

        mamba remove rgi


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

- `Python 3.6+ <https://www.python.org/>`_
- `NCBI BLAST 2.16.0 <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
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
- `samtools 1.21 <https://github.com/samtools/samtools>`_
- `bamtools 2.5.2 <https://github.com/pezmaster31/bamtools>`_
- `bedtools 2.31.1 <https://github.com/arq5x/bedtools2>`_
- `Jellyfish 2.2.10 <https://github.com/gmarcais/Jellyfish>`_
- `Bowtie2 2.5.4 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
- `BWA 0.7.18 <https://github.com/lh3/bwa>`_
- `KMA 1.4.17 <https://bitbucket.org/genomicepidemiology/kma/src/master>`_


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

