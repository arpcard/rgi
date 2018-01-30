Resistance Gene Identifier (RGI) 
--------------------------------------------

This application is used to predict resistome(s) from protein or nucleotide data based on homology and SNP models. The application uses data from `CARD database <https://card.mcmaster.ca/>`_.

Table of Contents
-------------------------------------

- `License`_
- `Requirements`_
- `Install Dependecies`_
- `Install RGI from project root`_
- `Running RGI Tests`_
- `Running RGI`_
- `Usage`_
- `Support & Bug Reports`_

License
--------
Use or reproduction of these materials, in whole or in part, by any non-academic organization whether or not for non-commercial (including research) or commercial purposes is prohibited, except with written permission of McMaster University. Commercial uses are offered only pursuant to a written license and user fee. To obtain permission and begin the licensing process, see `CARD website <https://card.mcmaster.ca/about>`_

Requirements
--------------------

- `Prodigal v2.6.3 <https://github.com/hyattpd/prodigal/wiki/Installation>`_
- `NCBI BLAST v2.6.0+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`_
- `DIAMOND v0.8.36 <https://ab.inf.uni-tuebingen.de/software/diamond>`_
- `Python 3.6 <https://www.python.org/>`_

Install Dependecies
--------------------

- pip3 install numpy
- pip3 install biopython
- pip3 install profilehooks
- pip3 install filetype
- pip3 install psycopg2
- pip3 install pandas

Install RGI from project root
-----------------------------

- pip3 install .


Running RGI Tests
-------------------

- cd tests
- pytest -v

Running RGI
-------------------

- rgi -h

Usage
-------------------

.. code-block:: sh

   usage: rgi <command> [<args>] 
            commands are:
               main    Runs rgi application
               tab     Creates a Tab-delimited from rgi results
               parser  Creates categorical .json files RGI wheel visualization. An input .json file containing the RGI results must be input.
               load    Loads CARD database json file
               clean   Removes BLAST databases and temporary files
               galaxy  Galaxy project wrapper

   Resistance Gene Identifier - 4.0.0

   positional arguments:
   command     Subcommand to run

   optional arguments:
   -h, --help  show this help message and exit

   Use the Resistance Gene Identifier to predict resistome(s) from protein or
   nucleotide data based on homology and SNP models. Check
   https://card.mcmaster.ca/download for software and data updates. Receive email
   notification of monthly CARD updates via the CARD Mailing List
   (https://mailman.mcmaster.ca/mailman/listinfo/card-l)


Support & Bug Reports
----------------------

Please log an issue on `github issue <https://github.com/arpcard/oop/issues>`_.

You can email the CARD curators or developers directly at `card@mcmaster.ca <mailto:card@mcmaster.ca>`_, via Twitter at `@arpcard <http://www.twitter.com/arpcard>`_.

