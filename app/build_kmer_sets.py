import os, sys, json, csv, argparse, glob
import app.make_kmer_json
from Bio import SeqIO, Seq
from argparse import RawTextHelpFormatter
from app.settings import APP_NAME, SOFTWARE_VERSION
"""
This script builds the CARD*k-mer sets.
Please provide location to the CARD*Resistomes&Variants nucleotide FASTAs
    and the location to the index file (index-for-model-sequences.txt).
"""

working_directory = os.path.join(os.getcwd())

def combine_variant_sequences(f1, f2, f3, f4):
    os.system(
        "cat {fasta_one} {fasta_two} {fasta_three} {fasta_four} > {prev_fasta}"
        .format(fasta_one=f1, fasta_two=f2,
        fasta_three=f3, fasta_four=f4,
        prev_fasta=os.path.join(working_directory, 'nucleotide_prevalence_all.fasta'))
    )

def split_variant_sequences(index, fasta):
    # Plasmids/Data type
    print('Getting CARD*Resistomes&Variants sequences metadata...')
    id_plas = set()
    id_chr = set()
    id_contig = set()
    id_path = {}
    with open(index, 'r') as index:
        i = csv.reader(index, delimiter='\t')
        for line in i:
            # Genomic
            if line[7] == 'ncbi_chromosome':
                if line[0] not in id_chr:
                    id_chr.add(line[0])
            elif line[7] == 'ncbi_plasmid':
                if line[0] not in id_plas:
                    id_plas.add(line[0])
            elif line[7] == 'ncbi_contig':
                if line[0] not in id_contig:
                    id_contig.add(line[0])

            # Taxonomy
            if line[0] not in id_path:
                id_path[line[0]] = [line[5]]
            elif line[0] in id_path:
                if line[5] not in id_path[line[0]]:
                    id_path[line[0]].append(line[5])
            else:
                print('Error')

    # Discard shared sequences in genomic sets
    p_and_chr = id_plas.intersection(id_chr)
    tplas = id_plas.difference(id_chr) # plasmid sequqences not in chromosomes
    tchr = id_chr.difference(id_plas)
    tcontig = id_contig.difference(tchr) # contig sequences not in chr
    tcontig = tcontig.difference(tplas)
    tcontig = tcontig.difference(p_and_chr)
    print('Writing genomic outputs...')
    with open(os.path.join(working_directory, 'both.fasta'), 'w') as both:
        with open(os.path.join(working_directory, 'plasmid.fasta'), 'w') as plasmids:
            with open(os.path.join(working_directory, 'chr.fasta'), 'w') as chrs:
                with open(os.path.join(working_directory, 'contig.fasta'), 'w') as contigs:
                    for entry in SeqIO.parse(fasta, 'fasta'):
                        prev_id = entry.id.split(":")[1].split("|")[0]
                        if prev_id in p_and_chr:
                            SeqIO.write(entry, both, 'fasta')
                        elif prev_id in tchr:
                            SeqIO.write(entry, chrs, 'fasta')
                        elif prev_id in tplas:
                            SeqIO.write(entry, plasmids, 'fasta')
                        elif prev_id in tcontig:
                            SeqIO.write(entry, contigs, 'fasta')
                        else:
                            print('Error')
    print('# of sequences in chr + plasmid:', len(p_and_chr))
    print('# of sequences in just plasmids', len(tplas))
    print('# of sequences in just chromosomes:', len(tchr))
    print('# of sequences in just contigs:', len(tcontig))

    # Discard prev sequences in multiple pathogens
    print('Writing taxonomic outputs...')
    species_count = 0
    genus_count = 0
    multi_count = 0
    with open(os.path.join(working_directory, 'species.fasta'), 'w') as single_path:
        with open(os.path.join(working_directory, 'genus.fasta'), 'w') as single_genus:
            with open(os.path.join(working_directory, 'multi.fasta'), 'w') as multi_path:
                for entry in SeqIO.parse(fasta, 'fasta'):
                    prev_id = entry.id.split(":")[1].split("|")[0]
                    if len(id_path[prev_id]) == 1:
                        # single[prev_id] = entry.seq
                        SeqIO.write(entry, single_path, 'fasta')
                        species_count += 1
                    else:
                        genus = [p.split()[0] for p in id_path[prev_id]]
                        if len(set(genus)) == 1:
                            SeqIO.write(entry, single_genus, 'fasta')
                            genus_count += 1
                        else:
                            # multi[prev_id] = entry.seq
                            SeqIO.write(entry, multi_path, 'fasta')
                            multi_count += 1
    print('# of sequences in single species:', species_count)
    print('# of sequences in single genus:', genus_count)
    print('# of sequences in multiple genus:', multi_count)

def count_kmers(k, threads):
    # Species
    print("Counting single species kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "species.temp.jf"),
                fasta=os.path.join(working_directory, "species.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    species_kmers = "species_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, species_kmers),
                db_path=os.path.join(working_directory, "species.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    # Genus
    print("Counting single genus kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "genus.temp.jf"),
                fasta=os.path.join(working_directory, "genus.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    genus_kmers = "genus_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, genus_kmers),
                db_path=os.path.join(working_directory, "genus.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    # Multi genus
    print("Counting multi genus kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "multi.temp.jf"),
                fasta=os.path.join(working_directory, "multi.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    multi_kmers = "multi_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, multi_kmers),
                db_path=os.path.join(working_directory, "multi.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    # Chr + Plasmid kmers
    print("Counting chr + plasmid kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "both.temp.jf"),
                fasta=os.path.join(working_directory, "both.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    both_kmers = "both_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, both_kmers),
                db_path=os.path.join(working_directory, "both.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    # Plasmid kmers
    print("Counting plasmid kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "plasmid.temp.jf"),
                fasta=os.path.join(working_directory, "plasmid.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    plasmid_kmers = "plasmid_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, plasmid_kmers),
                db_path=os.path.join(working_directory, "plasmid.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    # Chr kmers
    print("Counting chromosome kmers...")
    # Counts kmers
    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "chr.temp.jf"),
                fasta=os.path.join(working_directory, "chr.fasta"), threads=threads)
    )
    # Writes jellyfish output to text file
    chr_kmers = "chr_{k}.out".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, chr_kmers),
                db_path=os.path.join(working_directory, "chr.temp.jf"))
    )
    # os.system(
    #     "rm mer_counts.jf"
    # )

    return species_kmers, genus_kmers, multi_kmers, both_kmers, plasmid_kmers, chr_kmers

def is_tool(name):
    import distutils.spawn
    if distutils.spawn.find_executable(name) is not None:
        return True
    else:
        return False

def main(args):
    # check if jellyfish is installed
    # if is_tool("jellyfish") is False:
    #     print("Missing dependency: jellyfish.\nPlease install from https://github.com/gmarcais/Jellyfish/releases")
    #     exit("Debug")

    prevalence_directory = args.input_directory
    card_fasta = args.card_fasta
    k = args.k
    batch_size = args.batch_size

    files = glob.glob(os.path.join(prevalence_directory,"*"))
    for f in files:
        if "index" in f:
            index = f
        elif "nucleotide_fasta_protein_homolog_model_variants.fasta" in f:
            f1 = f
        elif "nucleotide_fasta_protein_overexpression_model_variants.fasta" in f:
            f2 = f
        elif "nucleotide_fasta_protein_variant_model_variants.fasta" in f:
            f3 = f
        elif "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta" in f:
            f4 = f
    # print(f1, f2, f3, f4)

    if not args.skip:
        try:
            index, f1, f2, f3, f4
        except NameError:
            exit("ERROR: cannot locate prevalence sequences and/or index-for-model-sequences.txt in directory")
        else:
            print("-- CONCATENATING ALL CARD*R&V SEQUENCES --")
            combine_variant_sequences(f1, f2, f3, f4)
            print("Done \n")

            print("-- SPLITTING CARD*R&V SEQUENCES --")
            split_variant_sequences(index, os.path.join(
                working_directory, "nucleotide_prevalence_all.fasta"))
            print("DONE \n")
    else:
        pass

    print("-- COUNTING KMERS --")
    species_kmers, genus_kmers, multi_kmers, both_kmers, plasmid_kmers, chr_kmers = count_kmers(k,args.threads)
    print("DONE \n")

    # """DEBUG"""
    # species_kmers = "species_temp.out"
    # genus_kmers = "genus_temp.out"
    # multi_kmers = "multi_temp.out"
    # plasmid_kmers = "plasmid_temp.out"
    # chr_kmers = "chr_temp.out"
    # both_kmers = "both_temp.out"

    species_file = os.path.join(working_directory, species_kmers)
    genus_file = os.path.join(working_directory, genus_kmers)
    multi_file = os.path.join(working_directory, multi_kmers)
    both_file = os.path.join(working_directory, both_kmers)
    plasmid_file = os.path.join(working_directory, plasmid_kmers)
    chr_file = os.path.join(working_directory, chr_kmers)
    variant_sequences = os.path.join(working_directory, \
                        "nucleotide_prevalence_all.fasta")

    print("-- CREATING KMER SETS --")
    app.make_kmer_json.make_json(plasmid_file, chr_file, both_file, genus_file, \
    species_file, multi_file, variant_sequences, index, k, args.threads, batch_size)
    print("DONE \n")

    print("-- CREATING AMR {}-MER SET --".format(k))

    os.system(
        "cat {card} {variants} > card_and_prevalence.fasta"
        .format(card=args.card_fasta, variants=variant_sequences)
    )

    os.system(
        "jellyfish count --mer-len {k} --threads {threads} --size 2G --output {jf_out} {fasta}"
        .format(k=k, jf_out=os.path.join(working_directory, "all_amr.temp.jf"),
                fasta=os.path.join(working_directory, "card_and_prevalence.fasta"), threads=args.threads)
    )

    amr_kmers = "all_amr_{k}mers.txt".format(k=k)
    os.system(
        "jellyfish dump --column --tab --output {name} {db_path}"
        .format(name=os.path.join(working_directory, amr_kmers),
                db_path=os.path.join(working_directory, "all_amr.temp.jf"))
    )

    print("Finished creating CARD*kmers set.")

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi kmer_build",
        description="{} - {} - Kmer Build \n\nBuilds the kmer sets for CARD*kmers".format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input_directory', dest="input_directory",
        help="input directory of prevalence data")
    parser.add_argument('-c', '--card', dest="card_fasta", required=True,
        help="fasta file of CARD reference sequences. If missing, run 'rgi card_annotation' to generate.")
    parser.add_argument('-k', dest="k", required=True,
        help="k-mer size (e.g., 61)")
    parser.add_argument('--skip', dest="skip", action='store_true',
        help="skips the concatenation and splitting of the CARD*R*V sequences.")
    parser.add_argument('-n','--threads', dest="threads", type=int,
            default=1, help="number of threads (CPUs) to use (default={})".format(1))
    parser.add_argument('--batch_size', dest='batch_size', type=int, default=100000, help='number of kmers to query at a time using pyahocorasick--the greater the number the more memory usage (default=100,000)')
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
