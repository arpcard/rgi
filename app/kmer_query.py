import os, sys, json, csv, argparse, multiprocessing, math
from app.settings import *
from collections import OrderedDict
from Bio import SeqIO, Seq

class CARDkmers(object):
    """
    Queries sequences against CARD*kmers
    """

    def __init__(self, input, bwt, rgi, fasta, k, threads, output, local, debug):
        # from arguments
        self.input_file = input
        self.bwt = bwt
        self.rgi = rgi
        self.fasta = fasta
        self.k = int(k)
        self.output = output
        self.local_database = local
        self.threads = threads

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
        if self.rgi:
            self.input_json_file = input
        if self.fasta:
            self.input_fasta_file = input

        # outputs
        self.working_directory = os.path.join(os.getcwd())
        self.base_name = os.path.basename(self.input_file)
        self.fasta_file = os.path.join(self.working_directory, "{}.fasta".format(self.base_name))
        self.output_json_file = os.path.join(self.working_directory, "{o}_{k}mer_analysis.json".format(o=self.output, k=self.k))

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
        else:
            print('error')
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

        logger.info("# of sequences in FASTA: {}".format(num_seq_total))
        logger.info("# of sequences with hits: {}".format(len(o_total)))
        logger.info("# of sequences too short: {}".format(short_total))
        logger.info("output json file: {}".format(self.output_json_file))
        print('done')
