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

