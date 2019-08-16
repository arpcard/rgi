import csv, re, argparse, multiprocessing, math, json, ahocorasick, os, sqlite3
from Bio import Seq, SeqIO
from app.settings import logger
import itertools
"""
This scripts creates the JSON to hold all kmer sets.
"""

def split_list(l,n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def query_kmers(temp_l, t, fasta, o, batch_size):
    """Finds kmers in variant sequences"""
    ### The fasta being referenced is "nucleotide_prevalence_all.fasta"
    list_size = math.ceil(len(temp_l)/(math.ceil(len(temp_l)/batch_size)))
    temp_l = (split_list(temp_l, list_size))
    # print("Total # of indexes:", len(temp_l))

    f = {}
    r = {}
    for i, l in enumerate(temp_l):
        try:

            logger.info('PROCESS {t}.{i}: Adding kmers'.format(t=t, i=i))
            A = ahocorasick.Automaton()
            for k in l:
                A.add_word(k,k)
            logger.info('PROCESS {t}.{i}: Making automaton'.format(t=t, i=i))
            A.make_automaton()

            seq_num = 0
            for entry in SeqIO.parse(fasta, "fasta"):
                seq_num += 1
                if seq_num % 25000 == 0:
                    logger.info('PROCESS {t}.{i}: Done querying {n} sequences'.format(t=t, n=seq_num, i=i))


                prev_id = entry.id.split(":")[1].split("|")[0]
                pathogens = id_path[prev_id] # get a list of pathogens that prev sequence is found in
                for tup in A.iter(str(entry.seq)):
                    kmer = tup[1]
                    rev_kmer = str(Seq.Seq(kmer).reverse_complement())
                    if kmer not in f:
                        f[kmer] = []
                    if rev_kmer not in r:
                        r[rev_kmer] = []
                    for p in pathogens:
                        if p not in f[kmer]:
                            f[kmer].append(p)
                        if p not in r[rev_kmer]:
                            r[rev_kmer].append(p)
                
            logger.info('PROCESS {t}.{i}: Done'.format(t=t, i=i))


        except Exception as e:
            logger.error(e)

    o.put((t, f, r))

def get_genomic_kmers(plasmid_file, chr_file, both_file):
    with open(plasmid_file, 'r') as f1:
        logger.info('Getting plasmid kmers...')
        p = csv.reader(f1, delimiter='\t')
        pkmers = {row[0] for row in p}
    with open(chr_file, 'r') as f2:
        logger.info('Getting chromosome kmers...')
        np = csv.reader(f2, delimiter='\t')
        nkmers = {row[0] for row in np}
    with open(both_file, 'r') as f3:
        logger.info('Getting both (plasmid + chr) kmers...')
        b = csv.reader(f3, delimiter='\t')
        bkmers = {row[0] for row in b}

        # Discard itnersect
        both_separate = pkmers.intersection(nkmers)
        # Remove non specific plasmid kmers
        pf = pkmers.difference(nkmers)
        pf = pf.difference(bkmers)
        # Remove non specific chromosome kmers
        cf = nkmers.difference(pkmers)
        cf = cf.difference(bkmers)
        for kmer in both_separate:
            bkmers.add(kmer)

        logger.info('Reverse complementing kmers...')
        # Reverse complement non-plasmid kmers
        rev_chr = {str(Seq.Seq(ck).reverse_complement()) for ck in cf}
        rev_plas = {str(Seq.Seq(pk).reverse_complement()) for pk in pf}
        rev_both = {str(Seq.Seq(bk).reverse_complement()) for bk in bkmers}

        # print(rev_chr)
        for false_hit in pf.intersection(rev_chr):
            # F kmer in plasmid is a R kmer in chr
            pf.remove(false_hit)
            bkmers.add(false_hit)
        for false_hit in pf.intersection(rev_both):
            pf.remove(false_hit)
            bkmers.add(false_hit)
        for false_hit in cf.intersection(rev_plas):
            cf.remove(false_hit)
            bkmers.add(false_hit)
        for false_hit in cf.intersection(rev_both):
            cf.remove(false_hit)
            bkmers.add(false_hit)

        return list(pf), list(cf), list(bkmers) ## type: sets

def get_taxon_kmers(single_file, multi_file, variant_sequences, index_file, k, type, threads, batch_size):
    """Gets taxonomic kmer sets"""
    global id_path
    id_path = {}
    logger.info('Getting pathogen metadata...')
    with open(index_file, 'r') as index:
        i = csv.reader(index, delimiter='\t')
        for line in i:
            if line[0] not in id_path:
                id_path[line[0]] = [line[5]]
            elif line[0] in id_path:
                if line[5] not in id_path[line[0]]:
                    id_path[line[0]].append(line[5])
            else:
                logger.error('Error')

    with open(single_file, 'r') as f1:
        logger.info('Getting single species kmers...')
        single = csv.reader(f1, delimiter='\t')
        skmers = {row[0] for row in single}
    with open(multi_file, 'r') as f2:
        logger.info('Getting multi species kmers...')
        multi = csv.reader(f2, delimiter='\t')
        mkmers = {row[0] for row in multi}

    logger.info('Discarding shared single- and multi-species kmers...')
    sf = list(skmers.difference(mkmers)) # kmers in only the single set
    logger.info('Remaining unique kmers: {}'.format(len(sf)))
    list_size = math.ceil(len(sf)/threads) # For 'n' threads
    l = list(split_list(sf, list_size))

    logger.info('Querying kmers...')
    # Multiprocessing
    output = multiprocessing.Queue()
    processes = []
    for ind in range(len(l)):
        process = multiprocessing.Process(target=query_kmers, args=(l[ind], ind, variant_sequences, output, batch_size))
        process.start()
        processes.append(process)
        # print(processes)

    results = [output.get() for process in processes]
    for process in processes:
        process.join()

    logger.info('Concatenating all data from all processes')
    f = {}
    r = {}
    for i in range(len(results)):
        f = {**f, **results[i][1]}
        r = {**r, **results[i][2]}

    to_delete = []
    # Filters kmers that are not unique
    for k in f:
        v = f[k]
        if type == "species":
            if len(v) > 1:
                to_delete.append(k)
        elif type == "genus":
            genus = [p.split()[0] for p in v]
            if len(set(genus)) > 1:
                to_delete.append(k)
        else:
            logger.error("Error")

    logger.info('before: {}'.format(len(f)))
    logger.info('to delete: {}'.format(len(to_delete)))
    for k in to_delete:
        rev_kmer = Seq.Seq(k).reverse_complement()
        del f[k]
        del r[rev_kmer]
    logger.info('after: {}'.format(len(f)))

    rev_multi = {Seq.Seq(mk).reverse_complement() for mk in mkmers}
    rev_temp = {Seq.Seq(tk).reverse_complement() for tk in to_delete}

    shared = 0
    same = 0
    single = 0
    forward = 0
    reverse = 0
    to_also_delete = []

    for k in r:
        if k in f:
            shared += 1
            if r[k] == f[k]:
                same += 1
            else: # Shared kmers to different pathogens in diff orientations
                single += 1
                to_also_delete.append(k)
    for k in to_also_delete:
        del f[k]
        del r[k]

    ts = set(f.keys())
    # print(len(ts.intersection(rev_temp)))
    for x in ts.intersection(rev_multi):
        del f[x]
        forward += 1
    for x2 in ts.intersection(rev_temp):
        del f[x2]
        forward += 1

    return f

def make_json(plasmid_file, chr_file, both_file, genus_file, species_file, \
                multi_file, variant_sequences, index_file, k, threads, batch_size):

    p, c, b = get_genomic_kmers(plasmid_file, chr_file, both_file) 
    s = get_taxon_kmers(species_file, multi_file, variant_sequences, index_file, k, "species", threads, batch_size)
    g = get_taxon_kmers(genus_file, multi_file, variant_sequences, index_file, k, "genus", threads, batch_size)

    final = {"p": p, "c": c, "b": b, "s": s, "g": g}

    with open(os.path.join(os.path.join(os.getcwd()), '{k}_kmer_db.json'.format(k=k)), 'w') as out:
        json.dump(final, out)


def main(args):
    plasmid_file = args.plasmid
    chr_file = args.chromosome
    both_file = args.both
    species_file = args.species
    genus_file = args.genus
    multi_file = args.multi
    variant_sequences = args.variants
    k = int(args.k)
    index_file = args.index
    threads = args.threads
    batch_size = args.batch_size

    make_json(plasmid_file, chr_file, both_file, genus_file, species_file, \
    multi_file, variant_sequences, index_file, k, threads, batch_size)

def run():
    parser = argparse.ArgumentParser(
        description='Creates a kmer catalogue')
    parser.add_argument('-p', dest="plasmid", required=True,
        help="Plasmid sequence kmers output from Jellyfish")
    parser.add_argument('-c', dest="chromosome", required=True,
        help="Chromosome sequence kmers output from Jellyfish")
    parser.add_argument('-b', dest="both", required=True,
        help="Plasmid + chromosome sequence kmers output from Jellyfish")
    parser.add_argument('-s', dest="species", required=True,
        help="Single species kmer output from Jellyfish")
    parser.add_argument('-g', dest="genus", required=True,
        help="Single genus kmer output from Jellyfish")
    parser.add_argument('-m', dest="multi", required=True,
        help="Multi species kmer output from Jellyfish")
    parser.add_argument('-v', dest="variants", required=True,
        help="CARD*Resistomes&Variants sequences (FASTA)")
    parser.add_argument('-i', dest="index", required=True,
        help="CARD*Resitomes&Variants index (index-for-model-sequences.txt)")
    parser.add_argument('-k', dest="k", required=True,
        help="kmer length")
    parser.add_argument('-n','--threads', dest="threads", type=int,
            default=1, help="number of threads (CPUs) to use (default={})".format(1))
    parser.add_argument('--batch_size', dest='batch_size', type=int, default=100000, help='Number of kmers to query at a time using pyahocorasick--the greater the number the more memory usage (default=100,000)')        
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
