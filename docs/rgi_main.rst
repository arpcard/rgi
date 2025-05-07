Analyzing Genomes, Genome Assemblies, Metagenomic Contigs, or Proteomes
-----------------------------------------------------------------------

If DNA FASTA sequences are submitted, RGI first predicts complete open reading frames (ORFs) using `Prodigal <https://github.com/hyattpd/Prodigal>`_ (ignoring those less than 30 bp) and analyzes the predicted protein sequences. This includes a secondary correction by RGI if Prodigal undercalls the correct start codon to ensure complete AMR genes are predicted. However, if Prodigal fails to predict an AMR ORF, RGI will produce a false negative result.

Short contigs, small plasmids, low quality assemblies, or merged metagenomic reads should be analyzed using Prodigal's algorithms for low quality/coverage assemblies (i.e. contigs <20,000 bp) and inclusion of partial gene prediction. If the low sequence quality option is selected, RGI uses Prodigal anonymous mode for open reading frame prediction, supporting calls of partial AMR genes from short or low quality contigs.

If protein FASTA sequences are submitted, RGI skips ORF prediction and uses the protein sequences directly.

The RGI analyzes genome or proteome sequences under a **Perfect**, **Strict**, and **Loose** (a.k.a. Discovery) paradigm. The Perfect algorithm is most often applied to clinical surveillance as it detects perfect matches to the curated reference sequences in CARD. In contrast, the Strict algorithm detects previously unknown variants of known AMR genes, including secondary screen for key mutations, using detection models with CARD's curated similarity cut-offs to ensure the detected variant is likely a functional AMR gene. The Loose algorithm works outside of the detection model cut-offs to provide detection of new, emergent threats and more distant homologs of AMR genes, but will also catalog homologous sequences and spurious partial matches that may not have a role in AMR. Combined with phenotypic screening, the Loose algorithm allows researchers to hone in on new AMR genes.

Within the **Perfect**, **Strict**, and **Loose** paradigm, RGI currently supports CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_, `protein variant models <https://card.mcmaster.ca/ontology/40293>`_, `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_, and `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_:

* **Protein Homolog Models** (PHM) detect protein sequences based on their similarity to a curated reference sequence, using curated BLASTP bitscore cut-offs, for example `NDM-1 <https://card.mcmaster.ca/ontology/36728>`_. Protein Homolog Models apply to all genes that confer resistance through their presence in an organism, such as the presence of a beta-lactamase gene on a plasmid. PHMs include a reference sequence and a bitscore cut-off for detection using BLASTP. A Perfect RGI match is 100% identical to the reference protein sequence along its entire length, a Strict RGI match is not identical but the bit-score of the matched sequence is greater than the curated BLASTP bit-score cutoff, Loose RGI matches have a bit-score less than the curated BLASTP bit-score cut-off.
* **Protein Variant Models** (PVM) perform a similar search as Protein Homolog Models (PHM), i.e. detect protein sequences based on their similarity to a curated reference sequence, but secondarily screen query sequences for curated sets of mutations to differentiate them from antibiotic susceptible wild-type alleles, for example `Acinetobacter baumannii gyrA conferring resistance to fluoroquinolones <https://card.mcmaster.ca/ontology/40507>`_. PVMs are designed to detect AMR acquired via mutation of house-keeping genes or antibiotic targets. PVMs include a protein reference sequence (often from antibiotic susceptible wild-type alleles), a curated bit-score cut-off, and mapped resistance variants. Mapped resistance variants may include any or all of single point mutations, insertions, or deletions curated from the scientific literature. A Strict RGI match has a BLASTP bit-score above the curated BLASTP cutoff value and contains at least one curated mutation from amongst the mapped resistance variants, while a Loose RGI match has a bit-score less than the curated BLASTP bit-score cut-off but still contains at least one curated mutation from amongst the mapped resistance variants.
* **Protein Overexpression Models** (POM) are similar to Protein Variant Models (PVM) in that they include a protein reference sequence, a curated BLASTP bitscore cut-off, and mapped resistance variants. Whereas PVMs are designed to detect AMR acquired via mutation of house-keeping genes or antibiotic targets, reporting only those with curated mutations conferring AMR, POMs are restricted to regulatory proteins and report both wild-type sequences and/or sequences with mutations leading to overexpression of efflux complexes, for example `MexS <https://card.mcmaster.ca/ontology/37193>`_. The former lead to efflux of antibiotics at basal levels, while the latter can confer clinical resistance. POMs include a protein reference sequence (often from wild-type alleles), a curated bit-score cut-off, and mapped resistance variants. Mapped resistance variants may include any or all of single point mutations, insertions, or deletions curated from the scientific literature. A Perfect RGI match is 100% identical to the wild-type reference protein sequence along its entire length, a Strict RGI match has a BLASTP bit-score above the curated BLASTP cutoff value may or may not contain at least one curated mutation from amongst the mapped resistance variants, while a Loose RGI match has a bit-score less than the curated BLASTP bit-score cut-off may or may not contain at least one curated mutation from amongst the mapped resistance variants.
* **Ribosomal RNA (rRNA) Gene Variant Models** (RVM) are similar to Protein Variant Models (PVM), i.e. detect  sequences based on their similarity to a curated reference sequence and secondarily screen query sequences for curated sets of mutations to differentiate them from antibiotic susceptible wild-type alleles, except RVMs are designed to detect AMR acquired via mutation of genes encoding ribosomal RNAs (rRNA), for example `Campylobacter jejuni 23S rRNA with mutation conferring resistance to erythromycin <https://card.mcmaster.ca/ontology/42445>`_. RVMs include a rRNA reference sequence (often from antibiotic susceptible wild-type alleles), a curated bit-score cut-off, and mapped resistance variants. Mapped resistance variants may include any or all of single point mutations, insertions, or deletions curated from the scientific literature. A Strict RGI match has a BLASTN bit-score above the curated BLASTN cutoff value and contains at least one curated mutation from amongst the mapped resistance variants, while a Loose RGI match has a bit-score less than the curated BLASTN bit-score cut-off but still contains at least one curated mutation from amongst the mapped resistance variants.

**Example**: The `Acinetobacter baumannii gyrA conferring resistance to fluoroquinolones <https://card.mcmaster.ca/ontology/40507>`_ Protein Variant Model has a bitscore cut-off of 1500 to separate **Strict** & **Loose** hits based on their similarity to the curated antibiotic susceptible reference protein AJF82744.1, but RGI will only report an antibiotic resistant version of this gene if the query sequence has the G79C or S81L substitutions:

.. image:: /images/gyrA.jpg

All RGI results are organized via the `Antibiotic Resistance Ontology <https://card.mcmaster.ca/ontology/36006>`_ classification: AMR Gene Family, Drug Class, and Resistance Mechanism. JSON files created at the command line can be `Uploaded at the CARD Website <https://card.mcmaster.ca/analyze/rgi>`_ for visualization, for example the Mycobacterium tuberculosis H37Rv complete genome (GenBank AL123456):

.. image:: /images/rgiwheel.jpg

**Note**: Users have the option of using BLAST or `DIAMOND <https://github.com/bbuchfink/diamond>`_ for generation of local alignments and assessment of bitscores within RGI. The default is BLAST, but DIAMOND generates alignments faster than BLAST and the RGI developers routinely assess DIAMOND's performance to ensure it calculates equivalent bitscores as BLAST given RGI's Perfect / Strict / Loose paradigm is dependant upon hand curated bitscore cut-offs. As such, RGI may not support the latest version of DIAMOND.

 > `What are CARD detection models and how are bitscore cut-offs determined? <https://github.com/arpcard/rgi/issues/140>`_

**UPDATED RGI version 6.0.0 onward: In earlier versions of RGI, by default all Loose matches of 95% identity or better were automatically listed as Strict, regardless of alignment length. At that time, this behaviour could only be suppressed by using the --exclude_nudge parameter. This default behaviour and the --exclude_nudge parameter have been discontinued. Loose matches of 95% identity or better can now only be listed (i.e., nudged) as Strict matches, regardless of alignment length, by use of the new --include_nudge parameter. As such, these often spurious results are no longer included in default RGI main output.**

Curation at CARD is routinely ahead of RGI software development, so not all parameters or models curated in CARD will be annotated in sequences analyzed using RGI. For example, RGI does not currently support CARD's `protein knockout models <https://card.mcmaster.ca/ontology/40354>`_, `protein domain meta-models <https://card.mcmaster.ca/ontology/40326>`_, `gene cluster meta-models <https://card.mcmaster.ca/ontology/40298>`_, or `efflux pump system meta-models <https://card.mcmaster.ca/ontology/41112>`_. In addition, while CARD's `protein variant models <https://card.mcmaster.ca/ontology/40293>`_, `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_, and `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_ are current supported by RGI, mutation screening currently only supports annotation of resistance-conferring SNPs via the `single resistance variant <https://card.mcmaster.ca/ontology/36301>`_ parameter. For example, here is a snapshot from CARD 4.0.0 for `protein variant models <https://card.mcmaster.ca/ontology/40293>`_:

+----------------------------------------------------------+------------------------------------------------+---------------------+
|    Parameters Among 242 PVMs                             | Frequency                                      | Supported by RGI    |
+==========================================================+================================================+=====================+
|    single resistance variant                             | 2398                                           |yes                  |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    nonsense mutation - Ter                               | 269                                            |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    multiple resistance variants                          | 114                                            |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    deletion mutation from nucleotide sequence            | 96                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    insertion mutation from nucleotide sequence           | 67                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    single resistance variant - Var                       | 61                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    snp in promoter region                                | 46                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    frameshift mutation - fs                              | 27                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    co-dependent single resistance variant                | 26                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    deletion mutation from peptide sequence               | 22                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    insertion mutation from peptide sequence              | 10                                             |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    co-dependent insertion/deletion - fs                  | 8                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    co-dependent single resistance variant - fs           | 8                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    co-dependent nonsense SNP - Ter                       | 5                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    co-dependent single resistance variant - Ter          | 5                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    insertion mutation                                    | 5                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    insertion mutation from peptide sequence - dup        | 4                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    snp in promoter region - Var                          | 3                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    disruptive mutation in regulatory element             | 2                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+
|    frameshift mutation - Ter                             | 1                                              |no                   |
+----------------------------------------------------------+------------------------------------------------+---------------------+

Lastly, analyzing metagenomic assemblies or merged metagenomic reads using RGI main is a computationally intensive approach, since each merged read or contig FASTA set may contain partial ORFs, requiring RGI to perform large amounts of BLAST/DIAMOND analyses against CARD reference proteins. However, this approach does (1) allow analysis of metagenomic sequences in protein space, overcoming issues of high-stringency read mapping relative to nucleotide reference databases (see below), and (2) allow inclusion of `protein variant models <https://card.mcmaster.ca/ontology/40293>`_, `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_, and `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_ when annotating the resistome (as outlined below, RGI bwt's read mapping algorithms do not support models that require screening for mutations).

 > `What RGI settings are best for a Metagenome-Assembled Genome (MAG)? <https://github.com/arpcard/FAQ#rgi-faqs>`_

Using RGI main
--------------

.. code-block:: sh

   rgi main -h

.. code-block:: sh

					usage: rgi main [-h] -i INPUT_SEQUENCE -o OUTPUT_FILE [-t {contig,protein}]
					                [-a {DIAMOND,BLAST}] [-n THREADS] [--include_loose]
					                [--include_nudge] [--local] [--clean] [--keep] [--debug]
					                [--low_quality] [-d {wgs,plasmid,chromosome,NA}] [-v]
					                [-g {PRODIGAL,PYRODIGAL}] [--split_prodigal_jobs]

					Resistance Gene Identifier - 6.0.2 - Main

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
					                        (default=16)
					  --include_loose       include loose hits in addition to strict and perfect
					                        hits (default: False)
					  --include_nudge       include hits nudged from loose to strict hits
					                        (default: False)
					  --local               use local database (default: uses database in
					                        executable directory)
					  --clean               removes temporary files (default: False)
					  --keep                keeps Prodigal CDS when used with --clean (default:
					                        False)
					  --debug               debug mode (default: False)
					  --low_quality         use for short contigs to predict partial genes
					                        (default: False)
					  -d {wgs,plasmid,chromosome,NA}, --data {wgs,plasmid,chromosome,NA}
					                        specify a data-type (default = NA)
					  -v, --version         prints software version number
					  -g {PRODIGAL,PYRODIGAL}, --orf_finder {PRODIGAL,PYRODIGAL}
					                        specify ORF finding tool (default = PRODIGAL)
					  --split_prodigal_jobs
					                        run multiple prodigal jobs simultaneously for contigs
					                        in a fasta file (default: False)


Loading CARD Reference Data for RGI main
`````````````````````````````````````````

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

Running RGI main with Genome or Assembly DNA Sequences
```````````````````````````````````````````````````````

The default settings for RGI main will include Perfect or Strict predictions via BLAST against CARD reference sequences for ORFs predicted by Prodigal from submitted nucleotide sequences, applying any additional mutation screening depending upon the detection model type, e.g. CARD's `protein homolog models <https://card.mcmaster.ca/ontology/40292>`_, `protein variant models <https://card.mcmaster.ca/ontology/40293>`_, `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_, and `protein over-expression models <https://card.mcmaster.ca/ontology/41091>`_. Prodigal ORF predictions will include complete start-to-stop ORFs only (ignoring those less than 30 bp).

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file 
        --local --clean

For AMR gene discovery, this can be expanded to include all Loose matches:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta
        --output_file /path/to/output_file --local --clean --include_loose

Or alternatively, users can select to list Loose matches of 95% identity or better as Strict matches, regardless of alignment length:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta
        --output_file /path/to/output_file --local --clean --include_nudge

Short contigs, small plasmids, low quality assemblies, or merged metagenomic reads should be analyzed using Prodigal's algorithms for low quality/coverage assemblies (i.e. contigs <20,000 bp) and inclusion of partial gene prediction. If the low sequence quality option is selected, RGI uses Prodigal anonymous mode for open reading frame prediction, supporting calls of partial AMR genes from short or low quality contigs:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta
        --output_file /path/to/output_file --local --clean --low_quality

Arguments can be used in combination. For example, analysis of metagenomic assemblies can be a computationally intensive approach so users may wish to use the faster DIAMOND algorithms, but the data may include short contigs with partial ORFs so the --low_quality flag may also be desirable. Partial ORFs may not pass curated bitscore cut-offs or novel samples may contain divergent alleles, so nudging 95% identity Loose matches to Strict matches may aid resistome annotation, although we suggest manual sorting of results by % identity or HSP length:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta
        --output_file /path/to/output_file --local --clean -a DIAMOND --low_quality
        --include_nudge

This same analysis can be threaded over many processors if high-performance computing is available:

   .. code-block:: sh

      rgi main --input_sequence /path/to/nucleotide_input.fasta
        --output_file /path/to/output_file --local --clean -a DIAMOND --low_quality
        --include_nudge --num_threads 40 --split_prodigal_jobs

Running RGI main with Protein Sequences
```````````````````````````````````````

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

If protein FASTA sequences are submitted, RGI skips ORF prediction and uses the protein sequences directly (thus excluding the `rRNA mutation models <https://card.mcmaster.ca/ontology/40295>`_). The same parameter combinations as above can be used, e.g. RGI annotating protein sequencing using the defaults:

   .. code-block:: sh

      rgi main --input_sequence /path/to/protein_input.fasta
        --output_file /path/to/output_file --local --clean -t protein

As above, for AMR gene discovery this can be expanded to include all Loose matches:

   .. code-block:: sh

      rgi main --input_sequence /path/to/protein_input.fasta
        --output_file /path/to/output_file --local --clean --include_loose -t protein

Other parameters can be used alone or in combination as above.

Running RGI main using GNU Parallel
````````````````````````````````````

System wide and writing log files for each input file. Note: add code below to script.sh then run with `./script.sh /path/to/input_files`.

   .. code-block:: sh

      #!/bin/bash
      DIR=`find . -mindepth 1 -type d`
      for D in $DIR; do
            NAME=$(basename $D);
            parallel --no-notice --progress -j+0 'rgi main -i {} -o {.} -n 16 -a diamond --clean --debug > {.}.log 2>&1' ::: $NAME/*.{fa,fasta};
      done

RGI main Tab-Delimited Output Details
`````````````````````````````````````

+----------------------------------------------------------+------------------------------------------------+
|    Field                                                 | Contents                                       |
+==========================================================+================================================+
|    ORF_ID                                                | Open Reading Frame identifier (internal to RGI)|
+----------------------------------------------------------+------------------------------------------------+
|    Contig                                                | Source Sequence                                |
+----------------------------------------------------------+------------------------------------------------+
|    Start                                                 | Start co-ordinate of ORF                       |
+----------------------------------------------------------+------------------------------------------------+
|    Stop                                                  | End co-ordinate of ORF                         |
+----------------------------------------------------------+------------------------------------------------+
|    Orientation                                           | Strand of ORF                                  |
+----------------------------------------------------------+------------------------------------------------+
|    Cut_Off                                               | RGI Detection Paradigm (Perfect, Strict, Loose)|
+----------------------------------------------------------+------------------------------------------------+
|    Pass_Bitscore                                         | Strict detection model bitscore cut-off        |
+----------------------------------------------------------+------------------------------------------------+
|    Best_Hit_Bitscore                                     | Bitscore value of match to top hit in CARD     |
+----------------------------------------------------------+------------------------------------------------+
|    Best_Hit_ARO                                          | ARO term of top hit in CARD                    |
+----------------------------------------------------------+------------------------------------------------+
|    Best_Identities                                       | Percent identity of match to top hit in CARD   |
+----------------------------------------------------------+------------------------------------------------+
|    ARO                                                   | ARO accession of match to top hit in CARD      |
+----------------------------------------------------------+------------------------------------------------+
|    Model_type                                            | CARD detection model type                      |
+----------------------------------------------------------+------------------------------------------------+
|    SNPs_in_Best_Hit_ARO                                  | Mutations observed in the ARO term of top hit  |
|                                                          | in CARD (if applicable)                        |
+----------------------------------------------------------+------------------------------------------------+
|    Other_SNPs                                            | Mutations observed in ARO terms of other hits  |
|                                                          | indicated by model id (if applicable)          |
+----------------------------------------------------------+------------------------------------------------+
|    Drug Class                                            | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
|    Resistance Mechanism                                  | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
|    AMR Gene Family                                       | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
|    Predicted_DNA                                         | ORF predicted nucleotide sequence              |
+----------------------------------------------------------+------------------------------------------------+
|    Predicted_Protein                                     | ORF predicted protein sequence                 |
+----------------------------------------------------------+------------------------------------------------+
|    CARD_Protein_Sequence                                 | Protein sequence of top hit in CARD            |
+----------------------------------------------------------+------------------------------------------------+
|    Percentage Length of Reference Sequence               | (length of ORF protein /                       |
|                                                          | length of CARD reference protein)              |
+----------------------------------------------------------+------------------------------------------------+
|    ID                                                    | HSP identifier (internal to RGI)               |
+----------------------------------------------------------+------------------------------------------------+
|    Model_id                                              | CARD detection model id                        |
+----------------------------------------------------------+------------------------------------------------+
|    Nudged                                                | TRUE = Hit nudged from Loose to Strict         |
+----------------------------------------------------------+------------------------------------------------+
|    Note                                                  | Reason for nudge or other notes                |
+----------------------------------------------------------+------------------------------------------------+
|    Hit_Start                                             | Start co-ordinate for HSP in CARD reference    |
+----------------------------------------------------------+------------------------------------------------+
|    Hit_End                                               | End co-ordinate for HSP in CARD reference      |
+----------------------------------------------------------+------------------------------------------------+
|    Antibiotic                                            | ARO Categorization                             |
+----------------------------------------------------------+------------------------------------------------+
|    AST_Source                                            | Source of antibiotic susceptibility data       |
|                                                          | (mutation-associated resistance only)          |
+----------------------------------------------------------+------------------------------------------------+

AST_Source: **Curated-R**: mutation data hand curated from the scientific literature, evaluated as conferring resistance (R). **CRyPTIC**: mutation data acquired from the `CRyPTIC catalog <https://pubmed.ncbi.nlm.nih.gov/35944069/>`_, evaluated as resistant (R), susceptible (S), or undetermined (U). **ReSeqTB**: mutation data acquired from the `ReSeqTB catalog <https://pubmed.ncbi.nlm.nih.gov/30337678/>`_, evaluated as conferring resistance (Minimal, Moderate, High), not conferring resistance (None), or Indeterminate. **WHO**: mutation data acquired from the `WHO 2023 catalog <https://www.who.int/publications/i/item/9789240082410>`_, evaluated as resistant (R), susceptible (S), or undetermined (U).

Generating Heat Maps of RGI main Results
````````````````````````````````````````

.. code-block:: sh

   rgi heatmap -h

.. code-block:: sh

				usage: rgi heatmap [-h] -i INPUT
				                   [-cat {drug_class,resistance_mechanism,gene_family}] [-f]
				                   [-o OUTPUT] [-clus {samples,genes,both}]
				                   [-d {plain,fill,text}] [--debug]

				Resistance Gene Identifier - 6.0.2 - Heatmap

				Creates a heatmap when given multiple RGI results.

				optional arguments:
				  -h, --help            show this help message and exit
				  -i INPUT, --input INPUT
				                        Directory containing the RGI .json files (REQUIRED)
				  -cat {drug_class,resistance_mechanism,gene_family}, --category {drug_class,resistance_mechanism,gene_family}
				                        The option to organize resistance genes based on a category.
				  -f, --frequency       Represent samples based on resistance profile.
				  -o OUTPUT, --output OUTPUT
				                        Name for the output EPS and PNG files.
				                        The number of files run will automatically
				                        be appended to the end of the file name.(default=RGI_heatmap)
				  -clus {samples,genes,both}, --cluster {samples,genes,both}
				                        Option to use SciPy's hiearchical clustering algorithm to cluster rows (AMR genes) or columns (samples).
				  -d {plain,fill,text}, --display {plain,fill,text}
				                        Specify display options for categories (deafult=plain).
				  --debug               debug mode

.. image:: /images/heatmap.jpg

RGI heatmap produces EPS and PNG image files. An example where rows are organized by AMR Gene Family and columns clustered by similarity of resistome is shown above.

Generate a heat map from pre-compiled RGI main JSON files, samples and AMR genes organized alphabetically:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
                --output /path/to/output_file

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by AMR gene family:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
                --output /path/to/output_file -cat gene_family -clus samples

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by Drug Class:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
                --output /path/to/output_file -cat drug_class -clus samples

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome and AMR genes organized by distribution among samples:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
                --output /path/to/output_file -clus both

Generate a heat map from pre-compiled RGI main JSON files, samples clustered by similarity of resistome (with histogram used for abundance of identical resistomes) and AMR genes organized by distribution among samples:

      .. code-block:: sh

            rgi heatmap --input /path/to/rgi_results_json_files_directory/
                --output /path/to/output_file -clus both -f

