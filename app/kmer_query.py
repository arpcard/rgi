import os, sys, json, csv, argparse, multiprocessing, math, pysam
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
        self.working_directory = os.path.join(os.getcwd())
        self.base_name = os.path.basename(self.input_file)
        self.fasta_file = os.path.join(self.working_directory, "{}.fasta".format(self.base_name))
        self.output_json_file = os.path.join(self.working_directory, "{o}_{k}mer_analysis.json".format(o=self.output, k=self.k))
        if self.bwt:
            self.input_bam_file = input
            self.bam_directory = os.path.dirname(self.input_bam_file)
            self.output_allele_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis.allele.txt".format(o=self.output, k=self.k))
            self.output_gene_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis.gene.txt".format(o=self.output, k=self.k))
        if self.rgi:
            self.input_rgi_file = input
            self.output_rgi_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis_rgi_summary.txt".format(o=self.output, k=self.k))
        if self.fasta:
            self.input_fasta_file = input
            self.output_fasta_summary = os.path.join(self.working_directory, "{o}_{k}mer_analysis_fasta_summary.txt".format(o=self.output, k=self.k))

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
        logger.info("start loading CARD*kmer {}-mer db".format(self.k))
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
        orf_list = []
        with open(self.input_rgi_file) as j:
            rgi_data = json.load(j)
        try:
            del rgi_data["_metadata"]
        except:
            pass

        with open(self.fasta_file, "w") as fasta:
            try:
                for key,value in rgi_data.items():
                    if isinstance(value, dict):
                        orf_id = key
                        orf_list.append(orf_id)
                        contig_id = orf_id.split()[0]
                        hsp = max(value.keys(), key=(lambda key: value[key]["bit_score"]))
                        model = value[hsp]["model_name"].replace(" ","_")
                        type_hit = value[hsp]["type_match"]
                        dna = value[hsp]["orf_dna_sequence"]

                        fasta.write(">{orf}__{hsp}__{model}__{type_hit}\n{dna}\n"
                                    .format(orf=contig_id, hsp=hsp, model=model, type_hit=type_hit, dna=dna))

            except Exception as e:
                print(e)

        return orf_list

    def get_bwt_sequences(self):
        """
        filter out unmapped and supplementary alignments.
        store RNAME, SAM flag, and MAPQ in header
        """

        os.system("samtools index {input_bam}".format(input_bam=self.input_bam_file))
        aligner = ""
        bam_file = pysam.AlignmentFile(self.input_bam_file, "rb")
        header = bam_file.text.split("\n")

        for h in header:
            if "@PG" in h:
                aligner = h.split("\t")[1]

        if aligner == "ID:KMA":
            os.system("""samtools view -F 4 -F 2048 {bam} | while read line; do awk '{cmd}'; done > {out}"""
                        .format(bam=self.input_bam_file, cmd="""{print ">"$1"__"$4"__"$3"__"$5"\\n"$11}""", out=self.fasta_file))
        else:
            os.system("""samtools view -F 4 -F 2048 {bam} | while read line; do awk '{cmd}'; done > {out}"""
                        .format(bam=self.input_bam_file, cmd="""{print ">"$1"__"$3"__"$2"__"$5"\\n"$10}""", out=self.fasta_file))



    def get_bwt_alignment_data(self, header):
        """
        bit-wise flag reference: http://blog.nextgenetics.net/?e=18
        """
        qname, model, flag, mapq = header.split("__")

        if flag.isdigit() is False:
            logger.error("failed to parse BAM file: {}, please check which aligner software was used to produce the BAM file in the RGI BWT step.".format(self.input_bam_file))
            exit()
        else:
            if int(flag) & 64: #
                read = "{}/1".format(qname)
            elif int(flag) & 128:
                read = "{}/2".format(qname)
            else: # single end data
                read = qname
            return read, model, flag, mapq

    # def get_rgi_data(self, header):
    #     orf, hsp, model, type_hit = header.split("__")
    #     return orf, hsp, model, type_hit

    def populate_rgi_json(self, orf, hsp, model, type_hit, num_kmers, amr_c, tax, gen, o):
        o[orf] = {'ORF': orf, 'contig': orf.split()[0], 'HSP': hsp, 'ARO_model': model.replace("_", " "), \
        'type_hit': type_hit, '#_of_kmers_in_sequence': num_kmers, \
        '#_of_AMR_kmers': amr_c, 'taxonomic_info': tax, 'genomic_info': gen}

    def populate_bwt_json(self, read, model, num_kmers, amr_c, flag, mapq, tax, gen, o):
        o[read] = {'reference': model, '#_of_kmers_in_sequence': num_kmers, \
        '#_of_AMR_kmers': amr_c, 'SAM_flag': int(flag), \
        'MAPQ': int(mapq), 'taxonomic_info': tax, \
        'genomic_info': gen}

    def populate_fasta_json(self, read, num_kmers, amr_c, tax, gen, o):
        o[read] = {'#_of_kmers_in_sequence': num_kmers, \
        '#_of_AMR_kmers': amr_c, 'taxonomic_info': tax, \
        'genomic_info': gen}

    def split_fasta(self, file):
        logger.info("getting sequences")
        iterator = SeqIO.parse(file, "fasta")
        temp_iterator = SeqIO.parse(file, "fasta")
        ns = sum(1 for i in temp_iterator) # counts number of sequences in generator
        list_size = math.ceil(ns/self.threads) # maximizes list size for threads available
        split_sequences = list(self.chunk_list(iterator, list_size)) if list_size > 0 else [[]] # returns a list of lists
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
                contig_id, hsp, model, type_hit = entry.id.split("__")
                for x in self.orf_list:
                    if x.split()[0] == contig_id:
                        orf_id = x
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
                        self.populate_rgi_json(orf_id, hsp, model, type_hit, num_kmers, amr_c, tax, gen, o)
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

    def single_species_bwt(self, read, path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele):
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    self.get_taxon_data(allele, sscp_allele, path)
                else:
                    self.get_taxon_data(allele, ssc_allele, path)
            else: # no chr kmers (only look at plasmid and genomic islands)
                if read['genomic_info']['plasmid'] > self.min-1:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, sscp_allele, path)
                    else:
                        # plasmid kmers present - plasmid
                        self.get_taxon_data(allele, ssp_allele, path)
                else:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        self.get_taxon_data(allele, sscp_allele, path)
                    else: # everything falls below cutoff
                        self.get_taxon_data(allele, ssu_allele, path)

        else:
            # no genomic kmers, but single species kmers
            self.get_taxon_data(allele, ssu_allele, path)

    def single_species_rgi(self, read, path):
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    prediction = "{} (chromosome or plasmid)".format(path)
                else:
                    prediction = "{} (chromosome)".format(path)
            else: # no chr kmers (only look at plasmid and genomic islands)
                if read['genomic_info']['plasmid'] > self.min-1:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        prediction = "{} (chromosome or plasmid)".format(path)
                    else:
                        # plasmid kmers present - plasmid
                        prediction = "{} (plasmid)".format(path)
                else:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        prediction = "{} (chromosome or plasmid)".format(path)
                    else: # everything falls below cutoff
                        prediction = "{} (no genomic info)".format(path)

        else:
            # no genomic kmers, but single species kmers
            prediction = "{} (no genomic info)".format(path)
        return prediction

    def single_genus_bwt(self, read, genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele):
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    # chr + plasmid kmers present - GI
                    self.get_taxon_data(allele, sgcp_allele, genus)
                else:
                    # chr kmers present - genomic
                    self.get_taxon_data(allele, sgc_allele, genus)
            else:
                if read['genomic_info']['plasmid'] > self.min-1:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, sgcp_allele, genus)
                    else:
                        # plasmid kmers present - plasmid
                        self.get_taxon_data(allele, sgp_allele, genus)
                else:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        self.get_taxon_data(allele, sgcp_allele, genus)
                    else:
                        self.get_taxon_data(allele, sgu_allele, genus)
        else:
            # no genomic kmers, but single genus kmers
            self.get_taxon_data(allele, sgu_allele, genus)

    def single_genus_rgi(self, read, genus):
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    # chr + plasmid kmers present - GI
                    prediction = "{} (chromosome or plasmid)".format(genus)
                else:
                    # chr kmers present - genomic
                    prediction = "{} (chromosome)".format(genus)
            else:
                if read['genomic_info']['plasmid'] > self.min-1:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        prediction = "{} (chromosome or plasmid)".format(genus)
                    else:
                        # plasmid kmers present - plasmid
                        prediction = "{} (plasmid)".format(genus)
                else:
                    if read['genomic_info']['chr + plasmid'] > self.min-1:
                        # chr + plasmid kmers present - GI
                        prediction = "{} (chromosome or plasmid)".format(genus)
                    else:
                        prediction = "{} (no genomic info)".format(genus)
        else:
            # no genomic kmers, but single genus kmers
            prediction = "{} (no genomic info)".format(genus)
        return prediction

    def ambiguous_bwt(self, read, allele, cp_allele, a_allele, m_allele, c_allele, taxon):
        # taxon - if
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    # differnt species and genus, but chr + plasmid kmer
                    self.get_ambiguous_data(allele, cp_allele)
                else:
                    # different species and genus only chr kmer (not mobile)
                    self.get_ambiguous_data(allele, c_allele)
            elif read['genomic_info']['plasmid'] > self.min-1:
                if read['genomic_info']['chr + plasmid'] > self.min-1:
                    # different species and genus, but chr + plasmid
                    self.get_ambiguous_data(allele, cp_allele)
                else:
                    # different species and genus, but plasmid
                    self.get_ambiguous_data(allele, m_allele)
            elif read['genomic_info']['chr + plasmid'] > self.min-1:
                # different species and genus, but chr + plasmid
                self.get_ambiguous_data(allele, cp_allele)
            else:
                # different species and genus, not enough genomic info
                if taxon:
                    self.get_ambiguous_data(allele, a_allele)
                # else:
                #     print("not enough info to make classification")
        else:
            if taxon:
                self.get_ambiguous_data(allele, a_allele)
            # else:
            #     print("not enough info to make classification")

    def ambiguous_rgi(self, read, taxon):
        # taxon - if
        if not(all(v == 0 for v in read['genomic_info'].values())):
            if read['genomic_info']['chr'] > self.min-1:
                if read['genomic_info']['plasmid'] + read['genomic_info']['chr + plasmid'] > self.min-1:
                    prediction = "Unknown taxonomy (chromosome or plasmid)"
                else:
                    # different species and genus only chr kmer (not mobile)
                    prediction = "Unknown taxonomy (chromosome)"
            elif read['genomic_info']['plasmid'] > self.min-1:
                if read['genomic_info']['chr + plasmid'] > self.min-1:
                    # different species and genus, but chr + plasmid
                    prediction = "Unknown taxonomy (chromosome or plasmid)"
                else:
                    # different species and genus, but plasmid
                    prediction = "Unknown taxonomy (plasmid)"
            elif read['genomic_info']['chr + plasmid'] > self.min-1:
                # different species and genus, but chr + plasmid
                prediction = "Unknown taxonomy (chromosome or plasmid)"
            else:
                # different species and genus, not enough genomic info
                if taxon:
                    prediction = "Unknown taxonomy and genomic context"
                else:
                    prediction = "N/A"
        else:
            if taxon:
                prediction = "Unknown taxonomy and genomic context"
            else:
                prediction = "N/A"
        return prediction

    def organize_summary_data(self, i, d, class_type):
        class_list = ["chromosome", "chromosome or plasmid", "plasmid",
                        "no genomic information"]
        hits = 0
        ss = []
        str = ""
        summary_str = ""
        for path in d[i]:
            hits += int(d[i][path])
            ss.append((path, int(d[i][path])))
        sorted_ss = sorted(ss, key = lambda x: x[1], reverse=True)
        for tup in sorted_ss:
            str =  str + "{p}: {c}; ".format(p=tup[0],c=tup[1])
            summary_str = summary_str + "{p} ({gen}): {c}; ".format(
                p=tup[0], gen=class_list[class_type], c=tup[1])
        return str, hits, summary_str


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

    def make_summaries(self, d, ssc, sscp, ssp, ssu, sgc, sgcp, sgp, sgu, cd, m, cp, ad):
        summary = []
        for a in d:
            total_hits = 0

            if a in ssc:
                ssc_str, ssc_hits, ssc_pred_str = self.organize_summary_data(a, ssc, 0)
            else:
                ssc_str, ssc_pred_str = "", ""
                ssc_hits = 0
            if a in sscp:
                sscp_str, sscp_hits, sscp_pred_str =self.organize_summary_data(a, sscp, 1)
            else:
                sscp_str, sscp_pred_str = "", ""
                sscp_hits = 0
            if a in ssp:
                ssp_str, ssp_hits, ssp_pred_str =self.organize_summary_data(a, ssp, 2)
            else:
                ssp_str, ssp_pred_str = "", ""
                ssp_hits = 0
            if a in ssu:
                ssu_str, ssu_hits, ssu_pred_str =self.organize_summary_data(a, ssu, 3)
            else:
                ssu_str, ssu_pred_str = "", ""
                ssu_hits = 0
            if a in sgc:
                sgc_str, sgc_hits, sgc_pred_str =self.organize_summary_data(a, sgc, 0)
            else:
                sgc_str, sgc_pred_str = "", ""
                sgc_hits = 0
            if a in sgcp:
                sgcp_str, sgcp_hits, sgcp_pred_str =self.organize_summary_data(a, sgcp, 1)
            else:
                sgcp_str, sgcp_pred_str = "", ""
                sgcp_hits = 0
            if a in sgp:
                sgp_str, sgp_hits, sgp_pred_str=self.organize_summary_data(a, sgp, 2)
            else:
                sgp_str, sgp_pred_str = "", ""
                sgp_hits = 0
            if a in sgu:
                sgu_str, sgu_hits, sgu_pred_str =self.organize_summary_data(a, sgu, 3)
            else:
                sgu_str, sgu_pred_str = "", ""
                sgu_hits = 0
            if a in cd:
                cc = cd[a]
            else:
                cc = 0
            if a in m:
                mc = m[a]
            else:
                mc = 0
            if a in cp:
                cpc = cp[a]
            else:
                cpc = 0
            if a in ad:
                ambc = ad[a]
            else:
                ambc = 0
            total_hits = ssc_hits + sscp_hits + ssp_hits + ssu_hits + sgc_hits + sgcp_hits + sgp_hits + sgu_hits + mc + cpc + ambc + cc

            prediction_str = ssc_pred_str + sscp_pred_str + ssp_pred_str + ssu_pred_str + sgc_pred_str + sgcp_pred_str + sgp_pred_str + sgu_pred_str


            summary.append({
                "id" : a,
                "mapped_reads_with_kmer_hits": total_hits,
                "prediction": prediction_str,
                'ssc': ssc_str,
                'sscp': sscp_str,
                'ssp': ssp_str,
                'ssu': ssu_str,
                'sgc': sgc_str,
                'sgcp': sgcp_str,
                'sgp': sgp_str,
                'sgu': sgu_str,
                'c': cc,
                'm': mc,
                'cp': cpc,
                'a': ambc,
                })

        return summary


    def parse_kmer_json(self, type):
        with open(self.output_json_file) as f:
            j = json.load(f)

        # counters = 0
        if type == "bwt":
            all_alleles = []
            ssc_allele = {} # key: ref seq, value: pathogen
            ssp_allele = {} # plasmid
            sscp_allele = {} # genomic island
            ssu_allele = {} # no genomic data
            sgc_allele = {} # single genus
            sgp_allele = {} # single genus plasmid
            sgcp_allele = {} # single genus genomic island
            sgu_allele = {}
            c_allele = {}
            m_allele = {} # promiscuous plasmid
            cp_allele = {} # genomic island
            a_allele = {} # ambiguous
            u = 0 # unknown - no hits/info
        elif type == "rgi":
            rgi_summary = []
        elif type == "fasta":
            fasta_summary = []

        hits = 0

        for read in j:
            if j[read]['#_of_kmers_in_sequence'] > 0:
                if type == "bwt":
                    allele = j[read]['reference']
                    if allele not in all_alleles:
                        all_alleles.append(allele)
                elif type == "rgi":
                    orf_id = j[read]["ORF"]
                    contig = j[read]["contig"]
                    model = j[read]["ARO_model"]
                    type_hit = j[read]["type_hit"]
                # elif type == "fasta":
                if j[read]['#_of_kmers_in_sequence'] > 0:
                    # single species kmers
                    if len(j[read]['taxonomic_info']['species']) == 1:
                        path = set(j[read]['taxonomic_info']['species'].keys()).pop()
                        if not j[read]['taxonomic_info']['genus']:
                            if j[read]['taxonomic_info']['species'][path] > self.min-1:
                                hits += 1
                                if type == "bwt":
                                    self.single_species_bwt(j[read], path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                else:
                                    prediction = self.single_species_rgi(j[read], path)
                            else:
                                if type == "bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                else:
                                    prediction =self.ambiguous_rgi(j[read], False)
                        else:
                            if len(j[read]['taxonomic_info']['genus']) == 1:
                                genus = set(j[read]['taxonomic_info']['genus'].keys()).pop()
                                if genus == path.split()[0]:
                                    if j[read]['taxonomic_info']['species'][path] + j[read]['taxonomic_info']['genus'][genus] > self.min-1:
                                    # genus of single species kmers match genus kmers
                                        hits += 1
                                        if type == "bwt":
                                            self.single_species_bwt(j[read], path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                        else:
                                            prediction = self.single_species_rgi(j[read], path)
                                    else:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], False)
                                else:
                                    if j[read]['taxonomic_info']['genus'][genus] > self.min-1 and j[read]['taxonomic_info']['species'][path] > self.min-1:
                                        hits += 1
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], True)
                                    elif j[read]['taxonomic_info']['genus'][genus] < self.min and j[read]['taxonomic_info']['species'][path] > self.min-1:
                                        hits += 1
                                        if type == "bwt":
                                            self.single_species_bwt(j[read], path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                        else:
                                            prediction = self.single_species_rgi(j[read], path)
                                    elif j[read]['taxonomic_info']['genus'][genus] > self.min-1  and j[read]['taxonomic_info']['species'][path] < self.min:
                                        hits += 1
                                        if type == "bwt":
                                            self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                        else:
                                            prediction = self.single_genus_rgi(j[read], genus)
                                    else:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], False)

                            else: # multiple genera
                                # at least 2 genera have kmer counts > 10
                                max_genus = max(j[read]['taxonomic_info']['genus'].keys(), key=(lambda key: j[read]['taxonomic_info']['genus'][key]))
                                if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['genus'].values()) > 1:
                                    hits += 1
                                    if type =="bwt":
                                        self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                    else:
                                        prediction = self.ambiguous_rgi(j[read], True)
                                elif j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                    if max_genus == path.split()[0]:
                                        if j[read]['taxonomic_info']['species'][path] + j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                        # genus of single species kmers match genus kmers
                                            hits += 1
                                            if type == "bwt":
                                                self.single_species_bwt(j[read], path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                            else:
                                                prediction = self.single_species_rgi(j[read], path)
                                        else:
                                            # this should never be called; debug
                                            print("this should never be called")
                                            if type =="bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], False)
                                    else:
                                        if j[read]['taxonomic_info']['species'][path] > self.min-1:
                                            hits += 1
                                            if type == "bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], True)
                                        elif j[read]['taxonomic_info']['species'][path] < self.min:
                                            hits += 1
                                            if type == "bwt":
                                                self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                            else:
                                                prediction = self.single_genus_rgi(j[read], genus)
                                else:
                                    if max_genus == path.split()[0]:
                                        if j[read]['taxonomic_info']['species'][path] + j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                            hits += 1
                                            if type =="bwt":
                                                self.single_species_bwt(j[read], path, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                            else:
                                                prediction = self.single_species_rgi(j[read], path)
                                        else:
                                            if type =="bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], False)
                                    else:
                                        if type =="bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], False)

                    elif len(j[read]['taxonomic_info']['species']) > 1:
                        hits += 1
                        max_species = max(j[read]['taxonomic_info']['species'].keys(), key=(lambda key: j[read]['taxonomic_info']['species'][key]))
                        # multiple single species
                        if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['species'].values()) > 1:
                            passing_species = [path for path in j[read]['taxonomic_info']['species'] if j[read]['taxonomic_info']['species'][path] > self.min-1]
                            genera = [path.split()[0] for path in passing_species]
                            if len(set(genera)) == 1:
                                genus = set(genera).pop()
                                if not j[read]['taxonomic_info']['genus']:
                                    if type == "bwt":
                                        self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                    else:
                                        prediction = self.single_genus_rgi(j[read], genus)
                                else:
                                    if len(j[read]['taxonomic_info']['genus']) == 1:
                                        if genus == set(j[read]['taxonomic_info']['genus'].keys()).pop():
                                            if type == "bwt":
                                                self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                            else:
                                                prediction = self.single_genus_rgi(j[read], genus)
                                        else:
                                            if type =="bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], True)
                                    else:
                                        max_genus = max(j[read]['taxonomic_info']['genus'].keys(), key=(lambda key: j[read]['taxonomic_info']['genus'][key]))
                                        if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['genus'].values()) > 1:
                                            if type =="bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], True)
                                        elif j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                            if max_genus == genus:
                                                # genus of single species kmers match genus kmers
                                                if type == "bwt":
                                                    self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                                else:
                                                    prediction = self.single_genus_rgi(j[read], genus)
                                            else:
                                                # we already know that the # of kmers for the species is > threshold
                                                if type == "bwt":
                                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                                else:
                                                    prediction = self.ambiguous_rgi(j[read], True)

                                        else:
                                            # aren't enough genus kmers
                                            if type == "bwt":
                                                self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                            else:
                                                prediction = self.single_genus_rgi(j[read], genus)
                            else:
                                if type =="bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                else:
                                    prediction = self.ambiguous_rgi(j[read], True)

                        elif j[read]['taxonomic_info']['species'][max_species] > self.min - 1:
                            if not j[read]['taxonomic_info']['genus']:
                                if type == "bwt":
                                    self.single_species_bwt(j[read], max_species, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                else:
                                    prediction = self.single_species_rgi(j[read], max_species)
                            else:
                                if len(j[read]['taxonomic_info']['genus']) == 1:
                                    genus = set(j[read]['taxonomic_info']['genus'].keys()).pop()
                                    if genus == max_species.split()[0]:
                                        if type == "bwt":
                                            self.single_species_bwt(j[read], max_species, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                        else:
                                            prediction = self.single_species_rgi(j[read], max_species)
                                    else:
                                        if j[read]['taxonomic_info']['genus'][genus] > self.min-1:
                                            if type == "bwt":
                                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                            else:
                                                prediction = self.ambiguous_rgi(j[read], True)
                                        elif j[read]['taxonomic_info']['genus'][genus] < self.min:
                                            if type == "bwt":
                                                self.single_species_bwt(j[read], max_species, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                            else:
                                                prediction = self.single_species_rgi(j[read], max_species)
                                else:
                                    max_genus = max(j[read]['taxonomic_info']['genus'].keys(), key=(lambda key: j[read]['taxonomic_info']['genus'][key]))
                                    # at least 2 genera have kmer counts > 10
                                    if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['genus'].values()) > 1:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], True)
                                    else:
                                        if j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                            if max_genus == max_species.split()[0]:
                                                if type == "bwt":
                                                    self.single_species_bwt(j[read], max_species, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                                else:
                                                    prediction = self.single_species_rgi(j[read], max_species)
                                            else:
                                                if type == "bwt":
                                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                                else:
                                                    prediction = self.ambiguous_rgi(j[read], True)
                                        elif j[read]['taxonomic_info']['genus'][max_genus] < self.min:
                                            if type == "bwt":
                                                self.single_species_bwt(j[read], max_species, allele, ssc_allele, ssp_allele, sscp_allele, ssu_allele)
                                            else:
                                                prediction = self.single_species_rgi(j[read], max_species)
                        else:
                            if not j[read]['taxonomic_info']['genus']:
                                if type == "bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                else:
                                    prediction = self.ambiguous_rgi(j[read], False)

                            else:
                                if len(j[read]['taxonomic_info']['genus']) == 1:
                                    genus = set(j[read]['taxonomic_info']['genus'].keys()).pop()
                                    if j[read]['taxonomic_info']['genus'][genus] > self.min-1:
                                        if type == "bwt":
                                            self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                        else:
                                            prediction = self.single_genus_rgi(j[read], genus)
                                    else:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                        else:
                                            self.ambiguous_rgi(j[read], False)
                                elif len(j[read]['taxonomic_info']['genus']) > 1:
                                    max_genus = max(j[read]['taxonomic_info']['genus'].keys(), key=(lambda key: j[read]['taxonomic_info']['genus'][key]))
                                    if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['genus'].values()) > 1:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], True)
                                    elif j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                        if type == "bwt":
                                            self.single_genus_bwt(j[read], max_genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                        else:
                                            prediction = self.single_genus_rgi(j[read], max_genus)
                                    else:
                                        if type == "bwt":
                                            self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                        else:
                                            prediction = self.ambiguous_rgi(j[read], False)

                    elif not j[read]['taxonomic_info']['species']:
                        # no species info, single genus info
                        if len(j[read]['taxonomic_info']['genus']) == 1:
                            hits += 1
                            genus = set(j[read]['taxonomic_info']['genus'].keys()).pop()
                            if j[read]['taxonomic_info']['genus'][genus] > self.min-1:
                                if type == "bwt":
                                    self.single_genus_bwt(j[read], genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                else:
                                    prediction = self.single_genus_rgi(j[read], genus)
                            else:
                                if type == "bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                else:
                                    prediction = self.ambiguous_rgi(j[read], False)
                        elif len(j[read]['taxonomic_info']['genus']) > 1:
                            hits += 1
                            max_genus = max(j[read]['taxonomic_info']['genus'].keys(), key=(lambda key: j[read]['taxonomic_info']['genus'][key]))
                            if sum(int(c) > self.min-1 for c in j[read]['taxonomic_info']['genus'].values()) > 1:
                                if type == "bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, True)
                                else:
                                    prediction = self.ambiguous_rgi(j[read], True)
                            elif j[read]['taxonomic_info']['genus'][max_genus] > self.min-1:
                                if type == "bwt":
                                    self.single_genus_bwt(j[read], max_genus, allele, sgcp_allele, sgc_allele, sgp_allele, sgu_allele)
                                else:
                                    prediction = self.single_genus_rgi(j[read], max_genus)
                            else:
                                if type == "bwt":
                                    self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                                else:
                                    prediction = self.ambiguous_rgi(j[read], False)
                        elif not j[read]['taxonomic_info']['genus']:
                            hits += 1
                            if type == "bwt":
                                self.ambiguous_bwt(j[read], allele, cp_allele, a_allele, m_allele, c_allele, False)
                            else:
                                prediction = self.ambiguous_rgi(j[read], False)
                        else:
                            print('error --->',  read)
                    else:
                        print('error --->',  read)

            if type == "rgi" or type == "fasta":
                tax_str = ""
                gen_str = ""
                for k,v in j[read]['taxonomic_info'].items():
                    for k2, v2 in v.items():
                        tax_str = tax_str + "{p}: {c}; ".format(p=k2, c=v2)
                for k,v in j[read]['genomic_info'].items():
                    gen_str = gen_str + "{g}: {c}; ".format(g=k, c=v)

                if type == "rgi":
                    rgi_summary.append({
                        "orf": orf_id,
                        "contig": contig,
                        "type_hit": type_hit,
                        "model": model,
                        "prediction": prediction,
                        "taxonomic_kmers": tax_str,
                        "genomic_kmers": gen_str
                        })
                elif type == "fasta":
                    fasta_summary.append({
                        "header": read,
                        "kmers": j[read]["#_of_kmers_in_sequence"],
                        "amr": j[read]["#_of_AMR_kmers"],
                        "prediction": prediction,
                        "taxonomic_kmers": tax_str,
                        "genomic_kmers": gen_str

                    })


        if type == "bwt":
            return all_alleles, ssc_allele, ssp_allele, sscp_allele, ssu_allele, sgc_allele, sgp_allele, sgcp_allele, sgu_allele, c_allele, m_allele, cp_allele, a_allele
        elif type == "rgi":
            return rgi_summary
        elif type == "fasta":
            return fasta_summary
        else:
            logger.error("error")


    def make_bwt_summary(self, all_alleles, ssc_allele, ssp_allele, sscp_allele, ssu_allele, sgc_allele, sgp_allele, sgcp_allele, sgu_allele, c_allele, m_allele, cp_allele, a_allele):

        # parse allele data for genes
        all_genes = []
        ssc_gene = {}
        ssp_gene = {}
        sscp_gene = {}
        ssu_gene = {}
        sgc_gene = {}
        sgp_gene = {}
        sgcp_gene = {}
        sgu_gene = {}
        c_gene = {}
        m_gene = {}
        cp_gene = {}
        a_gene = {}

        for a in all_alleles:
            gene = (a.split('|')[2].split(':')[1]).replace("_"," ")
            # print(gene)
            if gene not in all_genes:
                all_genes.append(gene)

        self.parse_taxonomy_alleles_to_genes(ssc_allele, ssc_gene)
        self.parse_taxonomy_alleles_to_genes(ssp_allele, ssp_gene)
        self.parse_taxonomy_alleles_to_genes(sscp_allele, sscp_gene)
        self.parse_taxonomy_alleles_to_genes(ssu_allele, ssu_gene)
        self.parse_taxonomy_alleles_to_genes(sgc_allele, sgc_gene)
        self.parse_taxonomy_alleles_to_genes(sgp_allele, sgp_gene)
        self.parse_taxonomy_alleles_to_genes(sgcp_allele, sgcp_gene)
        self.parse_taxonomy_alleles_to_genes(sgu_allele, sgu_gene)
        self.parse_non_taxonomy_alleles_to_genes(c_allele, c_gene)
        self.parse_non_taxonomy_alleles_to_genes(m_allele, m_gene)
        self.parse_non_taxonomy_alleles_to_genes(cp_allele, cp_gene)
        self.parse_non_taxonomy_alleles_to_genes(a_allele, a_gene)

        # make summaries
        allele_summary = self.make_summaries(all_alleles, ssc_allele, sscp_allele, \
            ssp_allele, ssu_allele, sgc_allele, sgcp_allele, sgp_allele, sgu_allele, c_allele, m_allele, cp_allele, a_allele)

        gene_summary = self.make_summaries(all_genes, ssc_gene, sscp_gene, \
            ssp_gene, ssu_allele, sgc_gene, sgcp_gene, sgp_gene, sgu_gene, c_gene, m_gene, cp_gene, a_gene)

        # write summaries to txt files
        with open(self.output_allele_summary, "w") as allele_output:
            writer = csv.writer(allele_output, delimiter='\t')
            writer.writerow([
                    'Reference Sequence',
                    'Mapped reads with kmer DB hits',
                    'CARD*kmer Prediction',
                    'Single species (chromosome) reads',
                    'Single species (chromosome or plasmid) reads',
                    'Single species (plasmid) reads',
                    'Single species (no genomic info) reads',
                    'Single genus (chromosome) reads',
                    'Single genus (chromosome or plasmid) reads',
                    'Single genus (plasmid) reads',
                    'Single genus (no genomic info) reads',
                    'Promiscuous plasmid reads',
                    'Unknown taxonomy (chromosome) reads',
                    'Unknown taxonomy (chromosome or plasmid) reads',
                    'Unknown taxonomy (no genomic info) reads',
                        ])
            for r in allele_summary:
                writer.writerow([
                    r['id'],
                    r["mapped_reads_with_kmer_hits"],
                    r["prediction"],
                    r['ssc'],
                    r['sscp'],
                    r['ssp'],
                    r['ssu'],
                    r['sgc'],
                    r['sgcp'],
                    r['sgp'],
                    r['sgu'],
                    r['m'],
                    r['c'],
                    r['cp'],
                    r['a']
                ])

        with open(self.output_gene_summary, "w") as gene_output:
            writer = csv.writer(gene_output, delimiter='\t')
            writer.writerow([
                    'ARO term',
                    'Mapped reads with kmer DB hits',
                    'CARD*kmer Prediction',
                    'Single species (chromosome) reads',
                    'Single species (chromosome or plasmid) reads',
                    'Single species (plasmid) reads',
                    'Single species (no genomic info) reads',
                    'Single genus (chromosome) reads',
                    'Single genus (chromosome or plasmid) reads',
                    'Single genus (plasmid) reads',
                    'Single genus (no genomic info) reads',
                    'Promiscuous plasmid reads',
                    'Unknown taxonomy (chromosome) reads',
                    'Unknown taxonomy (chromosome or plasmid) reads',
                    'Unknown taxonomy (no genomic info) reads',
                        ])
            for r in gene_summary:
                writer.writerow([
                    r['id'],
                    r["mapped_reads_with_kmer_hits"],
                    r["prediction"],
                    r['ssc'],
                    r['sscp'],
                    r['ssp'],
                    r['ssu'],
                    r['sgc'],
                    r['sgcp'],
                    r['sgp'],
                    r['sgu'],
                    r['m'],
                    r['c'],
                    r['cp'],
                    r['a']
                ])

    def make_rgi_summary(self, summary):
        with open(self.output_rgi_summary, "w") as rgi_output:
            writer = csv.writer(rgi_output, delimiter="\t")
            writer.writerow([
                    "ORF_ID",
                    "Contig",
                    "Cut_Off",
                    "Best_Hit_ARO",
                    "CARD*kmer Prediction",
                    "Taxonomic kmers",
                    "Genomic kmers"
                ])
            for r in summary:
                writer.writerow([
                    r["orf"],
                    r["contig"],
                    r["type_hit"],
                    r["model"],
                    r["prediction"],
                    r["taxonomic_kmers"],
                    r["genomic_kmers"]
                ])

    def make_fasta_summary(self, summary):
        with open(self.output_fasta_summary, "w") as fasta_output:
            writer = csv.writer(fasta_output, delimiter="\t")
            writer.writerow([
                    "Sequence",
                    "Total # kmers",
                    "# of AMR kmers",
                    "CARD*kmer Prediction",
                    "Taxonomic kmers",
                    "Genomic kmers",
                ])
            for r in summary:
                writer.writerow([
                    r["header"],
                    r["kmers"],
                    r["amr"],
                    r["prediction"],
                    r["taxonomic_kmers"],
                    r["genomic_kmers"],
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
            logger.info("input rgi results file: {}".format(self.input_rgi_file))
            self.orf_list = self.get_rgi_sequences()
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

        # print("# of sequences queried: {}".format(num_seq_total))
        # print("# of sequences with hits: {}".format(len(o_total)))
        # print("# of sequences too short: {}".format(short_total))
        # print("output json file: {}".format(self.output_json_file))
        # print('done querying')

        # generate text summaries

        print('creating kmer summaries')
        if self.bwt:
            all_alleles, ssc_allele, ssp_allele, sscp_allele, ssu_allele, sgc_allele, \
                sgp_allele, sgcp_allele, sgu_allele, c_allele, m_allele, cp_allele, a_allele = self.parse_kmer_json("bwt")

            self.make_bwt_summary(all_alleles, ssc_allele, ssp_allele, sscp_allele, \
                ssu_allele, sgc_allele, sgp_allele, sgcp_allele, sgu_allele, c_allele, m_allele, cp_allele, a_allele)
        elif self.rgi:
            rgi_summary = self.parse_kmer_json("rgi")
            self.make_rgi_summary(rgi_summary)
        elif self.fasta:
            fasta_summary = self.parse_kmer_json("fasta")
            self.make_fasta_summary(fasta_summary)
        print("done creating kmer summaries")
