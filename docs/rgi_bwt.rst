.. contents:: Page Contents
   :local:
   
Analyzing Metagenomic Reads
---------------------------

RGI can align short DNA sequences in FASTQ format using `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ , `BWA <http://bio-bwa.sourceforge.net>`_ , or `KMA <https://bitbucket.org/genomicepidemiology/kma/src/master>`_ against CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_. The default and recommended read aligner is `KMA <https://bitbucket.org/genomicepidemiology/kma/src/master>`_ due to its documented `better performance for redundant databases <https://pubmed.ncbi.nlm.nih.gov/30157759/>`_ such as CARD. While CARD is not truly redundant, i.e. there are no identical reference sequences, CARD does reflect the `AMR alelle network problem <https://pubmed.ncbi.nlm.nih.gov/29335005/>`_ in that many sequences are very similar. For example, the nucleotide sequences of TEM-1 and TEM-2 are `99% similar with no alignment gaps </images/TEM-alignment.jpg>`_. A sample generating short reads from a legitimate TEM-1 gene may result in reads aligned among TEM-1, TEM-2, or other TEM beta-lactamases depending upon the alignment algorithm chosen. The `KMA publication <https://pubmed.ncbi.nlm.nih.gov/30157759/>`_ and our own simulations find KMA best resolves this issue:

.. image:: /images/simulation.jpg
The above illustrates simulated 90x short read coverage from seven antibiotic resistance gene nucleotide reference sequences in CARD (catB, OXA-1, AAC(6')-Ib, NDM-1, BRP(MBL), QnrB1, CTX-M-15), subsequently aligned with RGI bwt against CARD using Bowtie2 or KMA algorithms. Reads are aligned to a single reference gene using KMA but for Bowtie2 the same reads are aligned across a selection of similar reference sequences, with associated lower MAPQ scores. Note that KMA has limits in its ability to resolve very similar sequences, e.g. all simulated catB3 reads were all aligned to catI and all simulated AAC(6')-Ib reads were aligned to AAC(6')-Ib-cr. These simulated data are available at: https://github.com/raphenya/read-mapping-analysis.

**UPDATED RGI version 6.0.0 onward: In earlier versions of RGI, by default RGI bwt aligned reads to reference sequences from CARD's protein homolog models, protein variant models, rRNA mutation models, and protein over-expression models. However, as outlined above, the latter three model types require comparison to CARD's curated lists of mutations known to confer phenotypic antibiotic resistance to differentiate alleles conferring resistance from antibiotic susceptible alleles, e.g. a wild-type gyrase susceptible to fluoroquinolones. As such, earlier versions of RGI were over-reporting antibiotic resistance genes by not checking for these curated mutations. For example, while the KMA algorithm reports SNPs relative to reference, RGI was not screening these SNPs against CARD. Read alignments against the protein variant model, rRNA mutation model, and protein over-expression model reference sequences can now only be listed by use of the new --include_other_models parameter, but at this time these results still do not include comparison to CARD's curated lists of mutations. As such, these often spurious results are no longer included in default RGI bwt output. Support for mutation screening models will be added to future versions of RGI bwt.**

For RGI bwt, FASTQ sequences can be aligned to the 'canonical' curated CARD reference sequences associated with the Antibiotic Resistance Ontology (i.e. sequences available in GenBank with clear experimental evidence of elevated MIC in a peer-reviewed journal available in PubMED) or additionally to the *in silico* predicted allelic variants available in CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ data set. The latter is highly recommended for non-clinical samples as the allelic diversity for AMR genes is greatly unrepresented in the published literature, with a strong bias towards clinical antibiotic resistance genes and pathogens, hampering high-stringency read mapping for samples with divergent alleles. Inclusion of CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ allows read mapping to predicted allelic variants and AMR gene homologs for a wide variety of pathogens, incorporation of CARD's `Prevalence Data <https://card.mcmaster.ca/prevalence>`_ for easier interpretation of predicted AMR genes, and ultimately use of k-mer classifiers for prediction of pathogen-of-origin for FASTQ reads predicted to encode AMR genes (see below).

 > `What data is included in CARD? Can I add unpublished data? <https://github.com/arpcard/FAQ#card-faqs>`_

CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ and `Prevalence Data <https://card.mcmaster.ca/prevalence>`_ (nicknamed WildCARD) were generated using the RGI to analyze molecular sequence data available in `NCBI Genomes <https://www.ncbi.nlm.nih.gov/genome/>`_ for hundreds of pathogens of interest (see `Sampling Table <https://card.mcmaster.ca/prevalence>`_). For each of these pathogens, complete chromosome sequences, complete plasmid sequences, genomic island sequences, and whole genome shotgun (WGS) assemblies were analyzed individually by RGI. RGI results were then aggregated to calculate prevalence statistics for distribution of AMR genes among pathogens and plasmids, predicted resistomes, and to produce a catalog of predicted AMR alleles. These data were predicted under RGI's **Perfect** and **Strict** paradigms (see above), the former tracking perfect matches at the amino acid level to the curated reference sequences and mutations in the CARD, while the latter predicts previously unknown variants of known AMR genes, including secondary screen for key mutations. The reported results are entirely dependant upon the curated AMR detection models in CARD, the algorithms available in RGI, the pathogens sampled, and the sequence data available at NCBI at their time of generation. RGI bwt will indicate if the reference sequence for aligned reads is from the 'canonical' curated CARD reference sequences or from CARD's Resistomes & Variants, allowing users to know if the underlying reference is an *in silico* prediction or experimentally validated resistance gene.

**Note**: While CARD's Resistomes & Variants increases the allelic diversity of the reference data for non-clinical samples, it does so at the cost of inflating the allele network problem outlined above. Summarizing results at the level of AMR Gene Family may be more accurate than summarizing at the level of individual antibiotic resistance genes.

**Note**: As RGI bwt makes no assumptions about pre-processing of metagenomics data, we suggest prior quality/adaptor trimming of reads with `skewer <https://github.com/relipmoc/skewer>`_ and deduplication of reads using `dedupe.sh <https://sourceforge.net/projects/bbmap/>`_. If needed, down-sampling of FASTQ data can be performed using `seqtk <https://github.com/lh3/seqtk>`_. Thanks to Allison Guitor of McMaster University for these suggestions.

Using RGI bwt
-------------

.. code-block:: sh

   rgi bwt -h

.. code-block:: sh

				usage: rgi bwt [-h] -1 READ_ONE [-2 READ_TWO] [-a {kma,bowtie2,bwa}]
				               [-n THREADS] -o OUTPUT_FILE [--debug] [--clean] [--local]
				               [--include_wildcard] [--include_other_models] [--include_baits]
				               [--mapq MAPQ] [--mapped MAPPED] [--coverage COVERAGE]

				Resistance Gene Identifier - 6.0.2 - BWT

				Aligns metagenomic reads to CARD and wildCARD reference using kma, bowtie2 or bwa and provide reports.

				optional arguments:
				  -h, --help            show this help message and exit
				  -1 READ_ONE, --read_one READ_ONE
				                        raw read one (qc and trimmed)
				  -2 READ_TWO, --read_two READ_TWO
				                        raw read two (qc and trimmed)
				  -a {kma,bowtie2,bwa}, --aligner {kma,bowtie2,bwa}
				                        select read aligner (default=kma)
				  -n THREADS, --threads THREADS
				                        number of threads (CPUs) to use (default=16)
				  -o OUTPUT_FILE, --output_file OUTPUT_FILE
				                        name of output filename(s)
				  --debug               debug mode (default=False)
				  --clean               removes temporary files (default=False)
				  --local               use local database (default: uses database in executable directory)
				  --include_wildcard    include wildcard (default=False)
				  --include_other_models
				                        include protein variant, rRNA variant, knockout, and protein overexpression models (default=False)
				  --include_baits       include baits (default=False)
				  --mapq MAPQ           filter reads based on MAPQ score (default=False)
				  --mapped MAPPED       filter reads based on mapped reads (default=False)
				  --coverage COVERAGE   filter reads based on coverage of reference sequence

**Note**: The mapq, mapped, and coverage filters are planned features and do not yet work (but values are reported for manual filtering). Support for AMR bait capture methods (--include_baits) is forthcoming.

`BWA <http://bio-bwa.sourceforge.net>`_ usage within RGI bwt:

   .. code-block:: sh

      bwa mem -M -t {threads} {index_directory} {read_one} > {output_sam_file}

`Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ usage within RGI bwt:

   .. code-block:: sh

      bowtie2 --very-sensitive-local --threads {threads} -x {index_directory}
        -U {unpaired_reads} -S {output_sam_file}

`KMA <https://bitbucket.org/genomicepidemiology/kma/src/master/>`_ usage within RGI bwt (default):

   .. code-block:: sh

      kma -mem_mode -ex_mode -1t1 -vcf -int {read_one} -t {threads}
        -t_db {index_directory} -o {output_sam_file}.temp -sam

Running RGI bwt with FASTQ files - Restricted to Protein Homolog Models
````````````````````````````````````````````````````````````````````````

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

As outlined above, metagenomics analyses may additionally include CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ protein homolog model reference data if desired. If you wish to include these reference data, additionally download the Resistomes & Variants (a.ka. WildCARD) data:

   .. code-block:: sh

      wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
      mkdir -p wildcard
      tar -xjf wildcard_data.tar.bz2 -C wildcard
      gunzip wildcard/*.gz

Pre-process the WildCARD reference data for metagenomics reads (note that the filenames *wildcard_database_v3.0.2.fasta* and *card_database_v3.0.1.fasta* plus the parameter *version_number* depend on the version of CARD data downloaded, please adjust accordingly):

   .. code-block:: sh

      rgi wildcard_annotation -i wildcard --card_json /path/to/card.json
        -v version_number > wildcard_annotation.log 2>&1
      rgi load --wildcard_annotation wildcard_database_v3.0.2.fasta
        --card_json /path/to/card.json
        --wildcard_index /path/to/wildcard/index-for-model-sequences.txt
        --card_annotation card_database_v3.0.1.fasta --local

RGI will use FASTQ files as provided, be sure to include linker and quality trimming, plus sorting or any other needed pre-processing prior to using RGI (see suggestions above). **Note**: RGI bwt will assume unpaired reads unless the -2 flag is used. The examples below assume paired reads.

The default settings for RGI bwt will align reads using KMA against CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_, i.e. reference sequences that do not require SNP mapping to predict resistance. The default uses only 'canonical' curated CARD reference sequences associated with the Antibiotic Resistance Ontology (i.e. sequences available in GenBank with clear experimental evidence of elevated MIC in a peer-reviewed journal available in PubMED):

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local

The same analysis can be expanded to use multiple processors:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local -n 20

Although not recommended (see above), an alternate read aligner can be used:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local -n 20 -a bowtie2

RGI bwt can use an expanded reference set by aligning reads to both 'canonical' CARD **and** CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ `WildCARD` variants:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local -n 20 --include_wildcard

Running RGI bwt with FASTQ files - All Model Types
```````````````````````````````````````````````````

RGI bwt can also be used to align reads to CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_ **plus** `protein variant models <https://card.mcmaster.ca/ontology/40293>`_, `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_, and `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_. As outlined above, the latter three model types require comparison to CARD's curated lists of mutations known to confer phenotypic antibiotic resistance to differentiate alleles conferring resistance from antibiotic susceptible alleles, but RGI bwt as of yet does not perform this comparison. Use these results with caution.

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

Also pre-process these reference data for metagenomics reads (note that the filename *card_database_v3.0.1.fasta* depends on the version of CARD data downloaded, please adjust accordingly). Note the use of the *_all* version of reference files when loading reference data for all model types:

   .. code-block:: sh

      rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1
      rgi load -i /path/to/card.json
        --card_annotation_all_models card_database_v3.0.1_all.fasta --local

As outlined above, metagenomics analyses may additionally include CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ protein homolog model reference data if desired. If you wish to include these reference data, additionally download the Resistomes & Variants (a.ka. WildCARD) data:

   .. code-block:: sh

      wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
      mkdir -p wildcard
      tar -xjf wildcard_data.tar.bz2 -C wildcard
      gunzip wildcard/*.gz

Pre-process the WildCARD reference data for metagenomics reads (note that the filenames *wildcard_database_v3.0.2.fasta* and *card_database_v3.0.1.fasta* plus the paramater *version_number* depend on the version of CARD data downloaded, please adjust accordingly). Note the use of the *_all* version of reference files when loading reference data for all model types:

   .. code-block:: sh

      rgi wildcard_annotation -i wildcard --card_json /path/to/card.json
        -v version_number > wildcard_annotation.log 2>&1
      rgi load --card_json /path/to/card.json
        --wildcard_annotation_all_models wildcard_database_v3.0.2_all.fasta
        --wildcard_index /path/to/wildcard/index-for-model-sequences.txt
        --card_annotation_all_models card_database_v3.0.1_all.fasta
        --local

RGI will use FASTQ files as provided, be sure to include linker and quality trimming, plus sorting or any other needed pre-processing prior to using RGI (see suggestions above). **Note**: RGI bwt will assume unpaired reads unless the -2 flag is used. The examples below assume paired reads.

The default settings for RGI bwt will align reads using KMA:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local --include_other_models

The same analysis can be expanded to use multiple processors:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local --include_other_models -n 20

Although not recommended (see above), an alternate read aligner can be used:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local --include_other_models -n 20 -a bowtie2

RGI bwt can use an expanded reference set by aligning reads to both 'canonical' CARD **and** CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ `WildCARD` variants:

   .. code-block:: sh

      rgi bwt --read_one /path/to/fastq/R1.fastq.gz
        --read_two /path/to/fastq/R2.fastq.gz --output_file output_prefix
        --local --include_other_models -n 20 --include_wildcard

RGI bwt Tab-Delimited Output Details
````````````````````````````````````

RGI bwt aligns FASTQ reads to the AMR alleles used as reference sequences, with results provided for allele mapping and summarized at the AMR gene level (i.e. summing allele level results by gene). Five tab-delimited files are produced:

+----------------------------------------------------------+------------------------------------------------+
|    File                                                  | Contents                                       |
+==========================================================+================================================+
|    output_prefix.allele_mapping_data.txt                 | RGI bwt read mapping results at allele level   |
+----------------------------------------------------------+------------------------------------------------+
|    output_prefix.gene_mapping_data.txt                   | RGI bwt read mapping results at gene level     |
+----------------------------------------------------------+------------------------------------------------+
|    output_prefix.artifacts_mapping_stats.txt             | Statistics for read mapping artifacts          |
+----------------------------------------------------------+------------------------------------------------+
|    output_prefix.overall_mapping_stats.txt               | Statistics for overall read mapping results    |
+----------------------------------------------------------+------------------------------------------------+
|    output_prefix.reference_mapping_stats.txt             | Statistics for reference matches               |
+----------------------------------------------------------+------------------------------------------------+

RGI bwt read mapping results at allele level
````````````````````````````````````````````

+----------------------------------------------------------+---------------------------------------------------+
|    Field                                                 | Contents                                          |
+==========================================================+===================================================+
|    Reference Sequence                                    | Reference allele to which reads have been mapped  |
+----------------------------------------------------------+---------------------------------------------------+
|    ARO Term                                              | ARO Term                                          |
+----------------------------------------------------------+---------------------------------------------------+
|    ARO Accession                                         | ARO Accession                                     |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Model Type                                  | CARD detection model type                         |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference DB                                          | Reference allele is from either CARD or WildCARD  |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Allele Source                               | See below                                         |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed in Genome(s)          | Has this allele sequence been observed in a CARD  |
|                                                          | Prevalence genome sequence?                       |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed in Plasmid(s)         | Has this allele sequence been observed in a CARD  |
|                                                          | Prevalence plasmid sequence?                      |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed Pathogen(s)           | CARD Prevalence pathogens bearing this allele     |
|                                                          | sequence. If Reference DB is CARD, pathogen used  |
|                                                          | as the reference in the CARD detection model will |
|                                                          | be shown. Use k-mers to verify pathogen-of-origin.|
+----------------------------------------------------------+---------------------------------------------------+
|    Completely Mapped Reads                               | Number of reads mapped completely to allele       |
+----------------------------------------------------------+---------------------------------------------------+
|    Mapped Reads with Flanking Sequence                   | Number of reads mapped incompletely to allele     |
+----------------------------------------------------------+---------------------------------------------------+
|    All Mapped Reads                                      | Sum of previous two columns                       |
+----------------------------------------------------------+---------------------------------------------------+
|    Percent Coverage                                      | Percent of reference allele covered by reads      |
+----------------------------------------------------------+---------------------------------------------------+
|    Length Coverage (bp)                                  | Base pairs of reference allele covered by reads   |
+----------------------------------------------------------+---------------------------------------------------+
|    Average MAPQ (Completely Mapped Reads)                | Average MAPQ value                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Mate Pair Linkage                                     | For mate pair sequencing, if a sister read maps to|
|                                                          | a different AMR gene, this is listed              |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Length                                      | Length (bp) of reference allele                   |
+----------------------------------------------------------+---------------------------------------------------+
|    AMR Gene Family                                       | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Drug Class                                            | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistance Mechanism                                  | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Depth                                                 | Depth of coverage (reported only when using KMA)  |
+----------------------------------------------------------+---------------------------------------------------+
|    SNPs                                                  | Single nucleotide polymorphisms observed from     |
|                                                          | mapped reads (reported only when using KMA and    |
|                                                          | with depth of at least 5).                        |
|                                                          | Not screened against curated SNPs in CARD.        |
+----------------------------------------------------------+---------------------------------------------------+
|    Consensus Sequence DNA                                | Nucleotide Consensus Sequence using mapped reads  |
|                                                          | (reported only when using KMA and                 |
|                                                          | with depth of at least 5).                        |
+----------------------------------------------------------+---------------------------------------------------+
|    Consensus Sequence Protein                            | Protein Consensus Sequence translated from DNA    |
|                                                          | (reported only when using KMA and                 |
|                                                          | with depth of at least 5).                        |
+----------------------------------------------------------+---------------------------------------------------+

**Reference Allele Source:**

Entries with *CARD Curation* are aligned to a reference allele from a published, characterized AMR gene, i.e. 'canonical CARD', and thus encode a 100% match to the reference protein sequence. Otherwise, entries will be reported as *in silico* allele predictions based on either **Perfect** or **Strict** RGI matches in CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_, with percent identity to the CARD reference protein reported. Matches with low values should be used with caution, as CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ has predicted a low identity AMR homolog.

RGI bwt read mapping results at gene level
``````````````````````````````````````````

+----------------------------------------------------------+---------------------------------------------------+
|    Field                                                 | Contents                                          |
+==========================================================+===================================================+
|    ARO Term                                              | ARO Term                                          |
+----------------------------------------------------------+---------------------------------------------------+
|    ARO Accession                                         | ARO Accession                                     |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Model Type                                  | CARD detection model type                         |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference DB                                          | Reference allele(s) are from CARD and/or WildCARD |
+----------------------------------------------------------+---------------------------------------------------+
|    Alleles with Mapped Reads                             | # of alleles for this AMR gene with mapped reads  |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Allele(s) Identity to CARD Reference Protein| See below                                         |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed in Genome(s)          | Have these allele sequences been observed in a    |
|                                                          | CARD Prevalence genome sequence?                  |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed in Plasmid(s)         | Have these allele sequences been observed in a    |
|                                                          | CARD Prevalence plasmid sequence?                 |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistomes & Variants: Observed Pathogen(s)           | CARD Prevalence pathogens bearing this allele     |
|                                                          | sequence. If Reference DB is CARD, pathogen used  |
|                                                          | as the reference in the CARD detection model will |
|                                                          | be shown. Use k-mers to verify pathogen-of-origin.|
+----------------------------------------------------------+---------------------------------------------------+
|    Completely Mapped Reads                               | Number of reads mapped completely to these alleles|
+----------------------------------------------------------+---------------------------------------------------+
|    Mapped Reads with Flanking Sequence                   | Number of reads mapped incompletely to these      |
|                                                          | alleles                                           |
+----------------------------------------------------------+---------------------------------------------------+
|    All Mapped Reads                                      | Sum of previous two columns                       |
+----------------------------------------------------------+---------------------------------------------------+
|    Average Percent Coverage                              | Average % of reference allele(s) covered by reads |
+----------------------------------------------------------+---------------------------------------------------+
|    Average Length Coverage (bp)                          | Average bp of reference allele(s) covered by reads|
+----------------------------------------------------------+---------------------------------------------------+
|    Average MAPQ (Completely Mapped Reads)                | Statistics for reference matches                  |
+----------------------------------------------------------+---------------------------------------------------+
|    Number of Mapped Baits                                | not yet supported                                 |
+----------------------------------------------------------+---------------------------------------------------+
|    Number of Mapped Baits with Reads                     | not yet supported                                 |
+----------------------------------------------------------+---------------------------------------------------+
|    Average Number of reads per Bait                      | not yet supported                                 |
+----------------------------------------------------------+---------------------------------------------------+
|    Number of reads per Bait Coefficient of Variation (%) | not yet supported                                 |
+----------------------------------------------------------+---------------------------------------------------+
|    Number of reads mapping to baits and mapping to       | not yet supported                                 |
|    complete gene                                         |                                                   |
+----------------------------------------------------------+---------------------------------------------------+
|    Number of reads mapping to baits and mapping to       | not yet supported                                 |
|    complete gene (%)                                     |                                                   |
+----------------------------------------------------------+---------------------------------------------------+
|    Mate Pair Linkage (# reads)                           | For mate pair sequencing, if a sister read maps to|
|                                                          | a different AMR gene, this is listed (# reads     |
|                                                          | supporting linkage in parentheses)                |
+----------------------------------------------------------+---------------------------------------------------+
|    Reference Length                                      | Length (bp) of reference sequences                |
+----------------------------------------------------------+---------------------------------------------------+
|    AMR Gene Family                                       | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Drug Class                                            | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+
|    Resistance Mechanism                                  | ARO Categorization                                |
+----------------------------------------------------------+---------------------------------------------------+

**Reference Allele(s) Identity to CARD Reference Protein:**

Gives range of *Reference Allele Source* values reported in the RGI bwt read mapping results at allele level, indicating the range of percent identity at the amino acid level of the encoded proteins to the corresponding CARD reference sequence. Matches with low values should be used with caution, as CARD's `Resistomes & Variants <https://card.mcmaster.ca/genomes>`_ has predicted a low identity AMR homolog.

