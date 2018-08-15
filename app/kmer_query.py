import os, sys, json, csv, argparse, multiprocessing, math
from app.settings import *
from collections import OrderedDict
from Bio import SeqIO, Seq

class CARDkmers(object):
    """
    Queries sequences against CARD*kmers
    """

    def __init__(self, input, bwt, rgi, fasta, k, min, threads, output, local, debug):
        # from arguments
        self.input_file = input
        self.bwt = bwt
        self.rgi = rgi
        self.fasta = fasta
        self.k = int(k)
        self.output = output
        self.local_database = local
        self.threads = threads
        self.min = int(min)

        # database
        self.db = path
        self.data = data_path

        if self.local_database:
            self.db = LOCAL_DATABASE
            self.data = LOCAL_DATABASE

        self.kmer_db = os.path.join(self.data, "{}mer_database.json".format(k))
        self.amr_kmers = os.path.join(self.data, "amr_{}mer.txt".format(k))

        # files
        if self.bwt:
            self.input_bam_file = input
            self.bam_directory = os.path.dirname(self.input_bam_file)
            # use bam directory to search for .txt idx stats output
        if self.rgi:
            self.input_json_file = input
        if self.fasta:
            self.input_fasta_file = input

        # outputs
        self.working_directory = os.path.join(os.getcwd())
        self.base_name = os.path.basename(self.input_file)
        self.fasta_file = os.path.join(self.working_directory, "{}.fasta".format(self.base_name))
        self.output_json_file = os.path.join(self.working_directory, "{o}_{k}mer_analysis.json".format(o=self.output, k=self.k))
        self.output_allele_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis.allele.txt".format(o=self.output, k=self.k))
        self.output_gene_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis.gene.txt".format(o=self.output, k=self.k))

        self.debug = debug
        if self.debug:
            logger.setLevel(10)

    def __repr__(self):
    	"""Returns CARDkmers class full object."""
    	return "CARDkmers({}".format(self.__dict__)

    def check_databases_exist(self):
        # Check if required files loaded
        if not os.path.exists(self.kmer_db):
            logger.error("Missing {}.".format(self.kmer_db))
            exit()
        if not os.path.exists(self.amr_kmers):
            logger.error("Missing {}.".format(self.amr_kmers))
            exit()

    def load_kmers(self):
        with open(self.kmer_db) as db:
            j = json.load(db)
        logger.info("loaded CARD*kmer {}-mer db successfully".format(self.k))
        with open(self.amr_kmers) as amr:
            reader = csv.reader(amr, delimiter='\t')
            amr_kmers = {row[0] for row in reader}
        logger.info("loaded all AMR {}-mer set successfully".format(self.k))
        return j, amr_kmers

    def get_rgi_sequences(self):
        logger.info("getting sequences")
        with open(self.input_json_file) as j:
            rgi_data = json.load(j)
        try:
            del rgi_data["_metadata"]
        except:
            pass

        with open(self.fasta_file, "w") as fasta:
            try:
                for key,value in rgi_data.items():
                    if isinstance(value, dict):
                        contig_id = key.split()[0]
                        hsp = max(value.keys(), key=(lambda key: value[key]["bit_score"]))
                        model = value[hsp]["model_name"].replace(" ","_")
                        type_hit = value[hsp]["type_match"]
                        dna = value[hsp]["orf_dna_sequence"]

                        fasta.write(">{node}__{hsp}__{model}__{type_hit}\n{dna}\n"
                                    .format(node=contig_id, hsp=hsp, model=model, type_hit=type_hit, dna=dna))

            except Exception as e:
                print(e)

    def get_bwt_sequences(self):
        """
        filter out unmapped and supplementary alignments.
        store RNAME, SAM flag, and MAPQ in header
        """
        os.system("""samtools view -F 4 -F 2048 {bam} | while read line; do awk '{cmd}'; done > {out}"""
                    .format(bam=self.input_bam_file, cmd="""{print ">"$1"__"$3"__"$2"__"$5"\\n"$10}""", out=self.fasta_file))

    def get_bwt_alignment_data(self, header):
        """
        bit-wise flag reference: http://blog.nextgenetics.net/?e=18
        """
        qname, model, flag, mapq = header.split("__")
        if int(flag) & 64: #
            read = "{}/1".format(qname)
        elif int(flag) & 128:
            read = "{}/2".format(qname)
        else: # single end data
            read = qname
        return read, model, flag, mapq

    def get_rgi_data(self, header):
        orf, hsp, model, type_hit = header.split("__")
        return orf, hsp, model, type_hit

    def populate_rgi_json(self, orf, hsp, model, type_hit, num_kmers, amr_c, tax, gen, o):
        o[orf] = {'ORF': orf, 'HSP': hsp, 'ARO_model': model.replace("_", " "), \
        'type_hit': type_hit, '# of kmers in sequence': num_kmers, \
        '# of AMR kmers': amr_c, 'taxonomic info': tax, 'genomic info': gen}

    def populate_bwt_json(self, read, model, num_kmers, amr_c, flag, mapq, tax, gen, o):
        o[read] = {'reference': model, '# of kmers in sequence': num_kmers, \
        '# of AMR kmers': amr_c, 'SAM flag': int(flag), \
        'MAPQ': int(mapq), 'taxonomic info': tax, \
        'genomic info': gen}

    def populate_fasta_json(self, read, num_kmers, amr_c, tax, gen, o):
        o[read] = {'# of kmers in sequence': num_kmers, \
        '# of AMR kmers': amr_c, 'taxonomic info': tax, \
        'genomic info': gen}

    def split_fasta(self, file):
        logger.info("getting sequences")
        iterator = SeqIO.parse(file, "fasta")
        temp_iterator = SeqIO.parse(file, "fasta")
        ns = sum(1 for i in temp_iterator) # counts number of sequences in generator
        list_size = math.ceil(ns/self.threads) # maximizes list size for threads available
        split_sequences = list(self.chunk_list(iterator, list_size)) # returns a list of lists
        logger.info("Using {t} threads to query {s} sequences".format(t=self.threads, s=ns))
        return split_sequences

    def chunk_list(self, iterator, n):
        """
        Reference: https://biopython.org/wiki/Split_large_file
        """
        # print(list(iterator))
        entry = True
        while entry:
            batch = []
            while len(batch) < n:
                try:
                    entry = next(iterator)
                except StopIteration:
                    entry = None
                if entry is None: # end of file
                    break
                batch.append(entry)
            if batch:
                yield batch


    def query_sequences(self, k, j, amr_kmers, fasta, type, *args):
        o = OrderedDict()
        ps = {k for k in j["p"]} # plamids
        bs = {k for k in j["b"]} # both
        cs = {k for k in j["c"]} # chr
        num_seq = 0
        short = 0

        for entry in fasta:
            num_seq += 1
            # if num_seq % 100000 == 0:
            #     logger.info('Done querying %d sequences' % num_seq)
            # o[entry.id] = {}
            amr_c = 0
            c = 0
            s = {} # species
            gd = {} # genus data
            pc = 0 # plasmids
            bc = 0 # both
            cc = 0 # chr
            dna = str(entry.seq)
            num_kmers = len(dna) - k + 1
            if type == 'rgi':
                orf, hsp, model, type_hit = self.get_rgi_data(entry.id)
            elif type == 'bwt':
                read, model, flag, mapq = self.get_bwt_alignment_data(entry.id)
            elif type == "fasta":
                read = entry.id
            else:
                logger.error("incompatible input type")
            if num_kmers > 1:
                for i in range(num_kmers):
                    kmer = dna[i:i+k]
                    if kmer in amr_kmers:
                        amr_c += 1
                    if kmer in j["s"]:
                        c += 1
                        if j["s"][kmer][0] not in s:
                            s[j["s"][kmer][0]] = 1
                        else:
                            s[j["s"][kmer][0]] += 1
                    elif kmer in j["g"]:
                        c += 1
                        gs = set([p.split()[0] for p in j["g"][kmer]])
                        # Check if only one genus in the list as intended
                        if len(gs) == 1:
                            genus = gs.pop()
                            if genus not in gd:
                                gd[genus] = 1
                            else:
                                gd[genus] += 1
                        elif len(gs) > 1:
                            print('Multiple genuses for this kmer', kmer, gs)
                        else:
                            print('Error, empty taxonomic data.')
                    if kmer in bs:
                        bc += 1
                    elif kmer in ps:
                        pc += 1
                    elif kmer in cs:
                        cc += 1
                # Populate json dictionary
                tax = OrderedDict({'species': {}, 'genus': {}})
                gen = OrderedDict({'chr + plasmid': bc, 'plasmid': pc, 'chr': cc})
                for species, sv in s.items():
                    tax['species'][species] = sv
                for genus,gv in gd.items():
                    tax['genus'][genus] = gv
                if all(not tax[k] for k in tax) and all(v == 0 for v in gen.values()):
                    pass
                else:
                    if type == "rgi":
                        self.populate_rgi_json(orf, hsp, model, type_hit, num_kmers, amr_c, tax, gen, o)
                    elif type == "bwt":
                        self.populate_bwt_json(read, model, num_kmers, amr_c, flag, mapq, tax, gen, o)
                    elif type == "fasta":
                        self.populate_fasta_json(read, num_kmers, amr_c, tax, gen, o)
                    else:
                        logger.error("incompatible input type")
            else:
                short += 1

        if type == "rgi":
            return num_seq, short, o
        elif type == "bwt" or type == "fasta":
            args[0].put((num_seq, short, o))
        else:
            logger.error("incompatible input type")

    def execute_threads(self, split_sequences, j, amr_kmers, type):
        output = multiprocessing.Queue()
        processes = []
        for ind in range(len(split_sequences)):
            process = multiprocessing.Process(target=self.query_sequences,
                args=(self.k, j, amr_kmers, split_sequences[ind], type, output))
            process.start()
            processes.append(process)
        results = [output.get() for process in processes]
        for process in processes:
            process.join()
        return results

    def get_taxon_data(self, name, d, tax): # name - gene or allele, d - dictionary, tax - species or genus name
        if name not in d:
            d[name] = {}
            if tax not in d[name]:
                d[name][tax] = 1
            else:
                d[name][tax] += 1
        else:
            if tax not in d[name]:
                d[name][tax] = 1
            else:
                d[name][tax] += 1

    def get_ambiguous_data(self, name, d):
        if name not in d:
            d[name] = 1
        else:
            d[name] += 1

    def single_species(self, read, path, allele, ss_allele, ssp_allele, ssgi_allele):
        if not(all(v == 0 for v in read['genomic info'].values())):
            if read['genomic info']['chr'] > 0:
                if read['genomic info']['plasmid'] + read['genomic info']['chr + plasmid'] > self.min-1:
                    self.get_taxon_data(allele, ssgi_allele, path)
                else:
                    self.get_taxon_data(allele, ss_allele, path)
            else: # no chr kmers (only look at plasmid and genomic islands)
                if read['genomic info']['plasmid'] > self.min-1:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, ssgi_allele, path)
                    else:
                        # plasmid kmers present - plasmid
                        self.get_taxon_data(allele, ssp_allele, path)
                else:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        self.get_taxon_data(allele, ssgi_allele, path)
                    else: # everything falls below cutoff
                        self.get_taxon_data(allele, ss_allele, path)

        else:
            # no genomic kmers, but single species kmers
            self.get_taxon_data(allele, ss_allele, path)

    def single_genus(self, read, genus, allele, sggi_allele, sg_allele, sgp_allele):
        if not(all(v == 0 for v in read['genomic info'].values())):
            if read['genomic info']['chr'] > 0:
                if read['genomic info']['plasmid'] + read['genomic info']['chr + plasmid'] > self.min-1:
                    # chr + plasmid kmers present - GI
                    self.get_taxon_data(allele, sggi_allele, genus)
                else:
                    # chr kmers present - genomic
                    self.get_taxon_data(allele, sg_allele, genus)
            else:
                if read['genomic info']['plasmid'] > self.min-1:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, sggi_allele, genus)
                    else:
                        # plasmid kmers present - plasmid
                        self.get_taxon_data(allele, sgp_allele, genus)
                else:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, sggi_allele, genus)
                    else:
                        self.get_taxon_data(allele, sg_allele, genus)
        else:
            # no genomic kmers, but single genus kmers
            self.get_taxon_data(allele, sg_allele, genus)

    def abiguous(self, read, allele, gi_allele, a_allele, m_allele):
        if not(all(v == 0 for v in read['genomic info'].values())):
            if read['genomic info']['chr'] > 0:
                if read['genomic info']['plasmid'] + read['genomic info']['chr + plasmid'] > self.min-1:
                    # differnt species and genus, but chr + plasmid kmer
                    self.get_ambiguous_data(allele, gi_allele)
                else:
                    # different species and genus only chr kmer (not mobile)
                    self.get_ambiguous_data(allele, a_allele)
            else:
                if read['genomic info']['plasmid'] > self.min-1:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        # different species and genus, but chr + plasmid
                        self.get_ambiguous_data(allele, gi_allele)
                    else:
                        # different species and genus, but plasmid
                        self.get_ambiguous_data(allele, m_allele)
                else:
                    if read['genomic info']['chr + plasmid'] > self.min-1:
                        # different species and genus, but chr + plasmid
                        self.get_ambiguous_data(allele, gi_allele)
                    else:
                        self.get_ambiguous_data(allele, a_allele)
        else:
            # different species and genus, no genomic info
            self.get_ambiguous_data(allele, a_allele)

    def organize_summary_data(self, i, d):
        hits = 0
        ss = []
        str = ""
        for path in d[i]:
            hits += int(d[i][path])
            ss.append((path, int(d[i][path])))
        sorted_ss = sorted(ss, key = lambda x: x[1], reverse=True)
        for tup in ss:
            str =  str + "{p}: {c};".format(p=tup[0],c=tup[1])
        return str, hits


    def parse_taxonomy_alleles_to_genes(self, ad, gd):
        for allele in ad:
            gene = (allele.split('|')[2].split(':')[1]).replace("_"," ")
            # print(gene)
            if gene not in gd:
                gd[gene] = ad[allele]
            else:
                for tax in ad[allele]:
                    if tax not in gd[gene]:
                        gd[gene][tax] = ad[allele][tax]
                    else:
                        gd[gene][tax] += ad[allele][tax]

    def parse_non_taxonomy_alleles_to_genes(self, ad, gd):
        for allele in ad:
            gene = allele.split('|')[2].split(':')[1].replace("_"," ")
            # print(gene)
            if gene not in gd:
                gd[gene] = ad[allele]
            else:
                gd[gene] += ad[allele]

    def make_summaries(self, d, ss, ssgi, ssp, sg, sggi, sgp, m, gi, ad):
        summary = []
        for a in d:
            total_hits = 0

            if a in ss:
                ss_str, ss_hits = self.organize_summary_data(a, ss)
            else:
                ss_str = ""
                ss_hits = 0
            if a in ssgi:
                ssgi_str, ssgi_hits =self.organize_summary_data(a, ssgi)
            else:
                ssgi_str = ""
                ssgi_hits = 0
            if a in ssp:
                ssp_str, ssp_hits =self.organize_summary_data(a, ssp)
            else:
                ssp_str = ""
                ssp_hits = 0
            if a in sg:
                sg_str, sg_hits =self.organize_summary_data(a, sg)
            else:
                sg_str = ""
                sg_hits = 0
            if a in sggi:
                sggi_str, sggi_hits =self.organize_summary_data(a, sggi)
            else:
                sggi_str = ""
                sggi_hits = 0
            if a in sgp:
                sgp_str, sgp_hits =self.organize_summary_data(a, sgp)
            else:
                sgp_str = ""
                sgp_hits = 0
            if a in m:
                mc = m[a]
            else:
                mc = 0
            if a in gi:
                gic = gi[a]
            else:
                gic = 0
            if a in ad:
                ambc = ad[a]
            else:
                ambc = 0
            total_hits = ss_hits + ssgi_hits + ssp_hits + sg_hits + sggi_hits + sgp_hits + mc + gic + ambc

            summary.append({
                "id" : a,
                "mapped_reads_with_kmer_hits": total_hits,
                'ss': ss_str,
                'ssgi': ssgi_str,
                'ssp': ssp_str,
                'sg': sg_str,
                'sggi': sggi_str,
                'sgp': sgp_str,
                'm': mc,
                'gi': gic,
                'a': ambc,
                })

        return summary

    def parse_bwt_kmer_json(self):
        with open(self.output_json_file) as f:
            j = json.load(f)

        # counters = 0
        all_alleles = []
        num_reads = 0
        ss_allele = {} # key: ref seq, value: pathogen
        ssp_allele = {} # plasmid
        ssgi_allele = {} # genomic island
        sg_allele = {} # single genus
        sgp_allele = {} # single genus plasmid
        sggi_allele = {} # single genus genomic island
        m_allele = {} # promiscuous plasmid
        gi_allele = {} # genomic island
        a_allele = {} # ambiguous
        u = 0 # unknown - no hits/info
        hits = 0

        for read in j:
            if j[read]['# of kmers in sequence'] > 0:
                num_reads += 1
                if j[read]['reference'] != '*': # not unmapped
                    allele = j[read]['reference']
                    if allele not in all_alleles:
                        all_alleles.append(allele)
                    # single species kmers
                    if len(j[read]['taxonomic info']['species']) == 1:
                        path = set(j[read]['taxonomic info']['species'].keys()).pop()
                        if not j[read]['taxonomic info']['genus']:
                            if j[read]['taxonomic info']['species'][path] > self.min-1:
                                hits += 1
                                self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                            else:
                                print("Not enough single species kmers to make the call.")

                        else:
                            if len(j[read]['taxonomic info']['genus']) == 1:
                                genus = set(j[read]['taxonomic info']['genus'].keys()).pop()
                                if genus == path.split()[0]:
                                    if j[read]['taxonomic info']['species'][path] + j[read]['taxonomic info']['genus'][genus] > self.min-1:
                                    # genus of single species kmers match genus kmers
                                        hits += 1
                                        self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                    else:
                                        print("Not enough single species or genus kmers")
                                else:
                                    if j[read]['taxonomic info']['genus'][genus] > self.min-1 and j[read]['taxonomic info']['species'][path] > self.min-1:
                                        hits += 1
                                        self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                                    elif j[read]['taxonomic info']['genus'][genus] < self.min and j[read]['taxonomic info']['species'][path] > self.min-1:
                                        hits += 1
                                        self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                    elif j[read]['taxonomic info']['genus'][genus] > self.min-1  and j[read]['taxonomic info']['species'][path] < self.min:
                                        hits += 1
                                        self.single_genus(j[read], genus, allele, sggi_allele, sg_allele, sgp_allele)
                                    else:
                                        print("Not enough single species or genus kmers")

                            else: # multiple genera
                                # at least 2 genera have kmer counts > 10
                                max_genus = max(j[read]['taxonomic info']['genus'].keys(), key=(lambda key: j[read]['taxonomic info']['genus'][key]))
                                if sum(int(c) > self.min-1 for c in j[read]['taxonomic info']['genus'].values()) > 1:
                                    hits += 1
                                    self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                                elif j[read]['taxonomic info']['genus'][max_genus] > self.min-1:
                                    if max_genus == path.split()[0]:
                                        if j[read]['taxonomic info']['species'][path] + j[read]['taxonomic info']['genus'][max_genus] > self.min-1:
                                        # genus of single species kmers match genus kmers
                                            hits += 1
                                            self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                        else:
                                            print("Not enough single species or genus kmers")
                                    else:
                                        if j[read]['taxonomic info']['species'][path] > self.min-1:
                                            hits += 1
                                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                                        elif j[read]['taxonomic info']['species'][path] < self.min:
                                            hits += 1
                                            self.single_genus(j[read], genus, allele, sggi_allele, sg_allele, sgp_allele)

                    elif len(j[read]['taxonomic info']['species']) > 1:
                        hits += 1
                        max_species = max(j[read]['taxonomic info']['species'].keys(), key=(lambda key: j[read]['taxonomic info']['species'][key]))
                        # multiple single species
                        if sum(int(c) > self.min-1 for c in j[read]['taxonomic info']['species'].values()) > 1:
                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                        elif j[read]['taxonomic info']['species'][max_species] > self.min - 1:
                            if not j[read]['taxonomic info']['genus']:
                                self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                            else:
                                if len(j[read]['taxonomic info']['genus']) == 1:
                                    genus = set(j[read]['taxonomic info']['genus'].keys()).pop()
                                    if genus == max_species.split()[0]:
                                        self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                    else:
                                        if j[read]['taxonomic info']['genus'][genus] > self.min-1:
                                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                                        elif j[read]['taxonomic info']['genus'][genus] < self.min:
                                            self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                else:
                                    max_genus = max(j[read]['taxonomic info']['genus'].keys(), key=(lambda key: j[read]['taxonomic info']['genus'][key]))
                                    # at least 2 genera have kmer counts > 10
                                    if sum(int(c) > self.min-1 for c in j[read]['taxonomic info']['genus'].values()) > 1:
                                        self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                                    elif j[read]['taxonomic info']['genus'][max_genus] > self.min-1:
                                        if max_genus == max_species.split()[0]:
                                            self.single_species(j[read], path, allele, ss_allele, ssp_allele, ssgi_allele)
                                        else:
                                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                        else:
                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)

                    elif not j[read]['taxonomic info']['species']:
                        # no species info, single genus info
                        if len(j[read]['taxonomic info']['genus']) == 1:
                            hits += 1
                            genus = set(j[read]['taxonomic info']['genus'].keys()).pop()
                            self.single_genus(j[read], genus, allele, sggi_allele, sg_allele, sgp_allele)

                        elif len(j[read]['taxonomic info']['genus']) > 1:
                            hits += 1
                            max_genus = max(j[read]['taxonomic info']['genus'].keys(), key=(lambda key: j[read]['taxonomic info']['genus'][key]))
                            if sum(int(c) > self.min-1 for c in j[read]['taxonomic info']['genus'].values()) > 1:
                                self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                            elif j[read]['taxonomic info']['genus'][max_genus] > self.min-1:
                                self.single_genus(j[read], genus, allele, sggi_allele, sg_allele, sgp_allele)
                            else:
                                self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)

                        elif not j[read]['taxonomic info']['genus']:
                            hits += 1
                            self.abiguous(j[read], allele, gi_allele, a_allele, m_allele)
                        else:
                            print('error --->',  read)
                    else:
                        print('error --->',  read)

        # parse allele data for genes
        all_genes = []
        ss_gene = {}
        ssp_gene = {}
        ssgi_gene = {}
        sg_gene = {}
        sgp_gene = {}
        sggi_gene = {}
        m_gene = {}
        gi_gene = {}
        a_gene = {}

        for a in all_alleles:
            gene = (a.split('|')[2].split(':')[1]).replace("_"," ")
            # print(gene)
            if gene not in all_genes:
                all_genes.append(gene)

        self.parse_taxonomy_alleles_to_genes(ss_allele, ss_gene)
        self.parse_taxonomy_alleles_to_genes(ssp_allele, ssp_gene)
        self.parse_taxonomy_alleles_to_genes(ssgi_allele, ssgi_gene)
        self.parse_taxonomy_alleles_to_genes(sg_allele, sg_gene)
        self.parse_taxonomy_alleles_to_genes(sgp_allele, sgp_gene)
        self.parse_taxonomy_alleles_to_genes(sggi_allele, sggi_gene)

        self.parse_non_taxonomy_alleles_to_genes(m_allele, m_gene)
        self.parse_non_taxonomy_alleles_to_genes(gi_allele, gi_gene)
        self.parse_non_taxonomy_alleles_to_genes(a_allele, a_gene)

        # make summaries
        allele_summary = self.make_summaries(all_alleles, ss_allele, ssgi_allele, \
            ssp_allele, sg_allele, sggi_allele, sgp_allele, m_allele, gi_allele, a_allele)
        gene_summary = self.make_summaries(all_genes, ss_gene, ssgi_gene, \
            ssp_gene, sg_gene, sggi_gene, sgp_gene, m_gene, gi_gene, a_gene)

        # write summaries to txt files
        with open(self.output_allele_summary, "w") as allele_output:
            writer = csv.writer(allele_output, delimiter='\t')
            writer.writerow([
                    'Reference Sequence',
                    'Mapped reads with kmer DB hits',
                    'Single species reads',
                    'Single species (GI) reads',
                    'Single species (P) reads',
                    'Single genus reads',
                    'Single genus (GI) reads',
                    'Single genus (P) reads',
                    'Promiscuous plasmid reads',
                    'Genomic island (GI) reads',
                    'Ambiguous reads',
                        ])
            for r in allele_summary:
                writer.writerow([
                    r['id'],
                    r["mapped_reads_with_kmer_hits"],
                    r['ss'],
                    r['ssgi'],
                    r['ssp'],
                    r['sg'],
                    r['sggi'],
                    r['sgp'],
                    r['m'],
                    r['gi'],
                    r['a']
                ])

        with open(self.output_gene_summary, "w") as gene_output:
            writer = csv.writer(gene_output, delimiter='\t')
            writer.writerow([
                    'ARO term',
                    'Mapped reads with kmer DB hits',
                    'Single species reads',
                    'Single species (GI) reads',
                    'Single species (P) reads',
                    'Single genus reads',
                    'Single genus (GI) reads',
                    'Single genus (P) reads',
                    'Promiscuous plasmid reads',
                    'Genomic island (GI) reads',
                    'Ambiguous reads',
                        ])
            for r in gene_summary:
                writer.writerow([
                    r['id'],
                    r["mapped_reads_with_kmer_hits"],
                    r['ss'],
                    r['ssgi'],
                    r['ssp'],
                    r['sg'],
                    r['sggi'],
                    r['sgp'],
                    r['m'],
                    r['gi'],
                    r['a']
                ])

    def run(self):
        # print args
        logger.info(json.dumps(self.__dict__, indent=2))

        logger.info("check for databases")
        self.check_databases_exist()

        # checks only one data type given
        if sum([self.bwt, self.rgi, self.fasta]) > 1:
            logger.error("Only specify one input type.")
            exit()

        # load kmers
        j, amr_kmers = self.load_kmers()

        if self.rgi:
            logger.info("input rgi results file: {}".format(self.input_json_file))
            self.get_rgi_sequences()
            iterator = SeqIO.parse(self.fasta_file, "fasta")
            num_seq_total, short_total, o_total = self.query_sequences(self.k, j, amr_kmers, iterator, "rgi")
            with open(self.output_json_file, "w") as oj:
                json.dump(o_total, oj)
        elif self.bwt:
            logger.info("input RGI*BWT bam file: {}".format(self.input_bam_file))
            self.get_bwt_sequences()
            # split sequences for threading
            split_sequences = self.split_fasta(self.fasta_file)
            # Threading
            results = self.execute_threads(split_sequences, j, amr_kmers, "bwt")

            o_total = {}
            num_seq_total = 0
            short_total = 0

            for i in range(len(results)):
                o_total.update(results[i][2])
                num_seq_total += results[i][0]
                short_total += results[i][1]

            with open(self.output_json_file, "w") as oj:
                json.dump(o_total, oj)
        elif self.fasta:
            logger.info("input fasta file: {}".format(self.input_fasta_file))
            # split sequences for threading
            split_sequences = self.split_fasta(self.input_fasta_file)
            # Threading
            results = self.execute_threads(split_sequences, j, amr_kmers, "fasta")

            o_total = {}
            num_seq_total = 0
            short_total = 0

            for i in range(len(results)):
                o_total.update(results[i][2])
                num_seq_total += results[i][0]
                short_total += results[i][1]

            with open(self.output_json_file, "w") as oj:
                json.dump(o_total, oj)
        else:
            logger.error("please specify an input file type")
            exit()

        print("# of sequences queried: {}".format(num_seq_total))
        print("# of sequences with hits: {}".format(len(o_total)))
        print("# of sequences too short: {}".format(short_total))
        print("output json file: {}".format(self.output_json_file))
        print('done querying')

        # generate text summaries
        self.parse_bwt_kmer_json()
        print("done creating allele and gene kmer summaries")
