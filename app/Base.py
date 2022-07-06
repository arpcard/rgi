#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""
Base classes used throughout RGI
"""
import os
from Bio import SeqIO
import json
from abc import ABCMeta, abstractmethod
from app.settings import logger
from Bio.Seq import Seq
from pyfaidx import Fasta

class RGIBase(object):
    """Interface for RGI"""
    __metaclass__ = ABCMeta

    @abstractmethod
    def from_string(self): pass

    @abstractmethod
    def from_args(self): pass

    @abstractmethod
    def validate_inputs(self): pass

    @abstractmethod
    def run(self): pass

    @abstractmethod
    def create_databases(self): pass

    @abstractmethod
    def run_blast(self): pass

    @abstractmethod
    def filter_process(self): pass

    @abstractmethod
    def output(self): pass


class BaseModel(object):
    def extract_nth_bar(self, bunch_str, n):
        """
        Parse the nth "|" delimited result field from aligner output
        and return the result information with the appropriate type

        Parameters
        ----------
        Args:
            bunch_str (str): "|" delimited result information from aligner
            n (int): "|" field index to extract

        Return:
            result (int, float, ascii byte-str): result with correct type
        """

        # to get appropriate field set offset
        start = n + 3
        end = n + 4

        # subset the string to find the relevant set of "|" delimited results
        subset_split_str = bunch_str.split('|')[start: end]

        temporary_str = "|".join(subset_split_str)

        # rebuild the string and extract the hit information after colon
        result = temporary_str[temporary_str.find(':')+2:]
        result = result.rstrip()

        # check if integer first
        if result.isdigit():
            return int(result)

        # otherwise try to coerce to float
        else:
            try:
                return float(result)

            # if that doesn't work encode into ascii-bytestring
            except ValueError:
                # return result.encode("ascii", "replace")
                return result

    def extract_nth_hash(self, bunch_str, n):
        """
        Parse and return the nth hash delimited field from
        alignment output

        Parameters
        ----------

        Args:
            bunch_str (str): '#' delimited string from aligner
            n (int): '#' field index to extract

        Return:
            result (int, str): extracted value from output in appropriate type
        """
        if "#" not in bunch_str:
            return 0
        else:
            arr = bunch_str.split("#")
            if n >= len(arr):
                return ""
            else:
                # if first two positional information fields return as integers
                if n == 1 or n == 2:
                    return int(arr[n])

                # strandedness data field
                elif n == 3:
                    # convert 1/0 strand indicator to +/- notation
                    if int(arr[n]) == 1:
                        return "+"
                    else:
                        return "-"

                # return as string if not specific strandedness or positional
                else:
                    return arr[n]

    def find_num_dash(self, subject, index):
        """
        Finds location of mutation by counting the
        number of dashes in the aligner subject output

        Parameters
        ----------

        Args:
            subject (str): aligner output string
            index (int): max output size

        Returns:
            dash_count (int): position of the mutation/SNP
        """
        dash_count = 0
        output = []

        for i in range(len(subject)):
            if subject[i] == '-':
                dash_count += 1
            else:
                output.append(subject[i])
            if len(output) == index:
                break

        return dash_count

    def get_submitted_protein_sequence(self, seq_filepath):
        """
        Parses sequence fasta into a dictionary keyed with the sequence IDs

        Parameters
        ----------

        Args:
            seq_filepath (str): sequence filepath

        Returns:
            submitted_proteins_dict (dict): dictionary of sequences of the
                                            format {seq_id: sequence string}
        """
        submitted_proteins_dict = {}

        if os.stat(seq_filepath).st_size != 0:
            for record in SeqIO.parse(seq_filepath, 'fasta'):
                submitted_proteins_dict[record.id] = str(record.seq)

        return submitted_proteins_dict

    def get_orf_dna_sequence(self, input_file, input_type):
        """
        Get the predicted open reading frame nucleotide sequence.

        Args:
            input_file (str): filepath of the input file
            input_type (str): [contig, read]
        """

        filename = os.path.basename(input_file)
        predicted_genes_dict = {}

        if input_type in ["contig"]:
            contig_filepath = os.path.join(self.working_directory,
                                           filename + ".temp.contigToORF.fsa")
            if os.stat(contig_filepath).st_size != 0:
                for record in SeqIO.parse(contig_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)

        elif input_type in ["read"]:
            read_filepath = os.path.join(self.working_directory,
                                         filename + ".temp.read.fsa")
            if os.stat(read_filepath).st_size != 0:
                for record in SeqIO.parse(read_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)
        else:
            raise ValueError("input_type invalid \
                    (must be 'contig' or 'read'): {}".format(input_type))

        # write json for all predicted file
        pjson = json.dumps(predicted_genes_dict)

        predicted_filepath = os.path.join(self.working_directory,
                                          filename + ".temp.predictedGenes.json")
        with open(predicted_filepath, 'w') as wf:
            wf.write(pjson)

        return predicted_genes_dict

    def get_orf_protein_sequence(self, input_file, input_type):
        """
        Get the predicted open reading frame nucleotide sequence.

        Args:
            input_file (str): filepath of the input file
            input_type (str): [contig, read]
        """

        filename = os.path.basename(input_file)
        predicted_genes_dict = {}

        if input_type in ["contig"]:
            contig_filepath = os.path.join(self.working_directory,
                                           filename + ".temp.contig.fsa")
            if os.stat(contig_filepath).st_size != 0:
                for record in SeqIO.parse(contig_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)

        elif input_type in ["read"]:
            read_filepath = os.path.join(self.working_directory,
                                         filename + ".temp.read.fsa")
            if os.stat(read_filepath).st_size != 0:
                for record in SeqIO.parse(read_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)
        else:
            raise ValueError("input_type invalid \
                    (must be 'contig' or 'read'): {}".format(input_type))

        # write json for all predicted file
        pjson = json.dumps(predicted_genes_dict)

        predicted_filepath = os.path.join(self.working_directory,
                                          filename + ".temp.predictedGenes.protein.json")
        with open(predicted_filepath, 'w') as wf:
            wf.write(pjson)

        return predicted_genes_dict

    def results(self, blast_results, query_id, perfect, strict , loose, include_nudge=False):
        """
        Sort results to perfect, strict, loose paradigm

        Parameters
        ----------

        Args:
            blast_results (dict): dictionary containing perfect, strict, and loose hits
            query_id (str): blast record query
            perfect (dict): dictionary containing perfect hits
            strict (dict): dictionary containing strict hits
            loose (dict): dictionary containing loose hits

        Returns:
            blast_results (dict): dictionary of sorted results
        """

        nudged = False
        if len(perfect) == 0 and len(strict) == 0 and len(loose) > 0:
            if include_nudge is True:
                nudged , loose = self.nudge_loose_to_strict(loose)

            if nudged is True and self.loose is False:
                blast_results[query_id] = loose
            elif self.loose is True:
                blast_results[query_id] = loose

        elif len(perfect) == 0:
            if len(strict) > 0:
                # TODO:: add try catch here
                try:
                    nudged , strict = self.nudge_strict_to_perfect(strict)
                    blast_results[query_id] = strict
                except Exception as e:
                    logger.error(e)
        else:
            if len(perfect) > 0:
                blast_results[query_id] = perfect

        return blast_results

    def nudge_strict_to_perfect(self, strict):
        """
        Nudge strict hits with missing n-terminus, c-terminus and alternate start codons

        Parameters
        ----------

        Args:
            strict (dict): dictionary containing strict hits

        Returns:
            nudged (bool): True or False
            strict (dict): dictionary containing strict or perfect hits
        """

        nudged = False

        # - check if there is 100% match with matching part to the reference
        # - getting matching protein then pull nucleotides from reference and translate
        # - check the start codons including alternates
        # - promote to perfect if the start codon is present in the N-terminus
        # 95 <= int(loose[i]["perc_identity"]) <= 100
        # int(strict[s]["perc_identity"]) == 100
        for s in strict:
            if int(strict[s]["perc_identity"]) == 100 and strict[s]["type_match"] == "Strict" and strict[s]["model_type_id"] in [40292] and self.input_type == 'contig':
                reference = strict[s]["sequence_from_broadstreet"]
                query = strict[s]["orf_prot_sequence"]
                orf_dna_sequence_ = ""
                orf_end_ = ""
                orf_start_ = ""
                _partial_protein = ""
                partial_protein = ""
                # Missing n-terminus or c-terminus
                if len(query) < len(reference) and query in reference:
                    length_nucleotides = (len(reference) - len(strict[s]["match"]))*3

                    # pull nucleotides from query or submitted sequence
                    partial_bases, orf_start_, orf_end_ = self.get_part_sequence(
                        self.input_sequence, strict[s]["orf_from"],
                        strict[s]["orf_start"], strict[s]["orf_end"],
                        length_nucleotides, strict[s]["orf_strand"],
                        strict[s]["ARO_name"]
                    )

                    # check if partial_bases are multiples of 3
                    if len(partial_bases) % 3 == 0:
                        # logger.debug("Missing part: [{}]".format(partial_bases))
                        if strict[s]["orf_strand"] == "-":
                            # update orf_end
                            strict[s]["orf_end_possible"] = orf_end_ + 1
                            strict[s]["orf_start_possible"] = strict[s]["orf_start"]
                            orf_dna_sequence_ = str(Seq(partial_bases).reverse_complement()) + strict[s]["orf_dna_sequence"]
                            partial_protein = str(Seq(partial_bases).reverse_complement().translate(table=11))
                            # logger.debug("Reverse strand: {}".format(partial_protein))
                        else:
                            # update orf_start
                            strict[s]["orf_start_possible"] = orf_start_
                            strict[s]["orf_end_possible"] = int(strict[s]["orf_end"]) + 1
                            orf_dna_sequence_ = partial_bases + strict[s]["orf_dna_sequence"]
                            partial_protein = str(Seq(partial_bases).translate(table=11)).strip("*")
                            # logger.debug("Forward strand: {}".format(partial_protein))

                        # logger.debug("Translated protein: [{}]".format(partial_protein))
                        # update start codon to M for all other alternate start codons
                        if len(partial_protein) > 0:
                            _partial_protein = partial_protein[0]
                            if partial_protein[0] in ["L","M","I","V"]:
                                _partial_protein = "M"+partial_protein[1:]

                        combine = _partial_protein + strict[s]["match"]

                        if combine == strict[s]["sequence_from_broadstreet"]:
                            nudged = True
                            # logger.debug("Missing n-terminus push to Perfect: {} #complete_gene {}".format(strict[s]["ARO_name"],
                            #     json.dumps({
                            #         "ARO_name": strict[s]["ARO_name"],
                            #         "strand": strict[s]["orf_strand"],
                            #         "header": strict[s]["orf_from"],
                            #         "orf_start_possible": strict[s]["orf_start_possible"],
                            #         "orf_end_possible": strict[s]["orf_end_possible"],
                            #         "orf_prot_sequence_possible": combine,
                            #         "orf_prot_sequence_possible_": partial_protein + strict[s]["match"],
                            #         "nudged": nudged
                            #     }, indent=2)
                            #     )
                            # )
                            strict[s]["type_match"] = "Perfect"
                            strict[s]["nudged"] = nudged
                            strict[s]["note"] = "possible complete gene, missing n-terminus"
                            strict[s]["missed_bases_by_prodigal"] = partial_bases
                            strict[s]["orf_dna_sequence_possible"] = orf_dna_sequence_
                            strict[s]["orf_prot_sequence_possible"] = combine

                    else:
                        pass
                        # logger.debug("incomplete nucleotides for gene: {} #partial_gene {}".format(strict[s]["ARO_name"],
                        #     json.dumps({
                        #                 "ARO_name": strict[s]["ARO_name"],
                        #                 "strand": strict[s]["orf_strand"],
                        #                 "header": strict[s]["orf_from"],
                        #                 "orf_start": strict[s]["orf_start"],
                        #                 "orf_end": strict[s]["orf_end"],
                        #                 "partial_bases": partial_bases
                        #             }, indent=2)
                        #     )
                        # )

                # reference contained within open reading frame
                elif len(query) > len(reference) and reference in query:
                    # check if the reference gene is complete in card db
                    if strict[s]["partial"] == "0":
                        # codons from table 11 https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG11
                        start_codon = strict[s]["dna_sequence_from_broadstreet"][:3]
                        stop_codon = strict[s]["dna_sequence_from_broadstreet"][-3:]
                        if start_codon in ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"] and stop_codon in ["TAA", "TAG", "TGA"]:
                            strict[s]["note"] = "possible complete gene, reference contained within open reading frame"
                            strict[s]["orf_prot_sequence_possible"] = ""

                            if strict[s]["orf_strand"] == "-":
                                strict[s]["orf_start_possible"] = strict[s]["orf_start"]
                                strict[s]["orf_end_possible"] = int(strict[s]["orf_start"]) + len(strict[s]["dna_sequence_from_broadstreet"]) - 1
                            else:
                                strict[s]["orf_start_possible"] = int(strict[s]["orf_end"]) - len(strict[s]["dna_sequence_from_broadstreet"]) + 1
                                strict[s]["orf_end_possible"] = int(strict[s]["orf_end"])

                            # pull nucleotides from query or submitted sequence
                            partial_bases, orf_start_, orf_end_ = self.get_part_sequence(
                                self.input_sequence, strict[s]["orf_from"],
                                strict[s]["orf_start_possible"], strict[s]["orf_end_possible"],
                                0, strict[s]["orf_strand"],
                                strict[s]["ARO_name"]
                            )

                            if strict[s]["orf_strand"] == "-":
                                strict[s]["orf_dna_sequence_possible"] = str(Seq(partial_bases).reverse_complement())
                            else:
                                strict[s]["orf_dna_sequence_possible"] = partial_bases

                            if len(strict[s]["orf_dna_sequence_possible"]) % 3 == 0:
                                orf_prot_sequence_possible = str(Seq(strict[s]["orf_dna_sequence_possible"]).translate(table=11)).strip("*")
                                strict[s]["orf_prot_sequence_possible"] = orf_prot_sequence_possible
                                if orf_prot_sequence_possible == strict[s]["sequence_from_broadstreet"]:
                                    strict[s]["type_match"] = "Perfect"
                                    nudged = True
                            else:
                                logger.warning("incorrect open reading frame for coordinate: {}-{} on strand {} for {}".format
                                (strict[s]["orf_start_possible"], strict[s]["orf_end_possible"], strict[s]["orf_strand"], strict[s]["orf_from"]))

                            strict[s]["nudged"] = nudged

                            # logger.debug("Reference contained within open reading frame push to Perfect: {} #complete_gene {}".format(strict[s]["ARO_name"],
                            #     json.dumps({
                            #         "ARO_name": strict[s]["ARO_name"],
                            #         "strand": strict[s]["orf_strand"],
                            #         "header": strict[s]["orf_from"],
                            #         "orf_start_possible": strict[s]["orf_start_possible"],
                            #         "orf_end_possible": strict[s]["orf_end_possible"],
                            #         "orf_dna_sequence_possible": strict[s]["orf_dna_sequence_possible"],
                            #         "orf_prot_sequence_possible": strict[s]["orf_prot_sequence_possible"],
                            #         "orf_start_": orf_start_,
                            #         "orf_end_": orf_end_,
                            #         "nudged": strict[s]["nudged"]
                            #     }, indent=2)
                            # ))
                        else:
                            logger.warning("partial gene in card reference: {}, start_codon: {},  stop_codon: {}".format(
                                strict[s]["ARO_name"],
                                start_codon,
                                stop_codon
                            ))
                    else:
                        # print(json.dumps(strict[s], indent=2))
                        logger.warning("partial gene in card reference: {}".format(
                            strict[s]["ARO_name"]
                        ))
                # orf and reference are overlapping
                '''
                elif reference not in query and query not in reference:
                    logger.warning("TODO:: orf and reference are overlapping")
                '''
            elif 95 <= int(strict[s]["perc_identity"]) < 100 and strict[s]["type_match"] == "Strict" and strict[s]["model_type_id"] in [40292] and self.input_type == 'contig':
                # logger.debug("95 <= 100, strict hit")
                # logger.debug(json.dumps(strict[s], indent=2))
                pass

        return nudged, strict

    def get_part_sequence(self, fasta_file, header, start, stop, nterminus, strand, name):
        """
        Pull part sequence from fasta file
        # https://github.com/mdshw5/pyfaidx
        # pip install pyfaidx

        Parameters
        ----------

        Args:
            fasta_file (str): input fasta file
            header (str): header for fasta sequence
            start (str): start coordinate
            stop (str): stop coordinate
            nterminus (int): length of missing sequence
            strand (str): strand
            name (str): gene name

        Returns:
            sequence (str): portion on a sequence
        """
        # remove the last 2 characters from header as this is appended by prodigal
        header = header[:header.rfind("_")]
        genes = False
        # logger.info("[PARTIAL] ARO: {} | contig: {} | filename: {}".format(name, header, fasta_file))
        try:
            genes = Fasta(fasta_file, sequence_always_upper=False, read_long_names=False, one_based_attributes=True)
        except Exception as e:
            logger.error(e)
        # logger.info(genes.records)
        if genes:
            # logger.debug(json.dumps({"strand":strand, "start":start, "stop":stop, "nterminus":nterminus}, indent=2))
            if strand == "-":
                _start = stop + 1
                _stop = stop + nterminus
                # logger.debug("grep sequence from {}|-|{}-{}".format(header,_start, _stop,))
                if nterminus == 0:
                    # logger.debug("grep sequence from {}|-|{}-{}".format(header,start, stop,))
                    return str(genes.get_spliced_seq( header, [[start, stop]])), start, stop
                else:
                    return str(genes.get_spliced_seq( header, [[_start, _stop]])), _start, _stop
            elif strand == "+":
                _start = start - nterminus
                _stop = start - 1
                if _start <= 0:
                    _start = 1
                if _stop <= 0:
                    _stop = 1
                # logger.debug("grep sequence from {}|+|{}-{}".format(header,_start, _stop))
                if nterminus == 0:
                    # logger.debug("grep sequence from {}|+|{}-{}".format(header,start, stop))
                    return str(genes.get_spliced_seq( header, [[start, stop]])), start, stop
                else:
                    return str(genes.get_spliced_seq( header, [[_start, _stop]])), _start, _stop

    def nudge_loose_to_strict(self, loose):
        """
        Nudge loose hits with at least 95 percent identity to be strict hits

        Parameters
        ----------

        Args:
            loose (dict): dictionary containing loose hits

        Returns:
            nudged (bool): True or False
            loose (dict): dictionary containing loose or strict hits
        """
        nudged = False
        # check if there are any loose hits that might be strict
        for i in loose:
            if 95 <= int(loose[i]["perc_identity"]) <= 100:
                # add to strict
                # logger.debug("loose hit with at least 95 percent identity push to Strict: {} #partial_gene".format(loose[i]["ARO_name"]))
                loose[i]["type_match"] = "Strict"
                loose[i]["nudged"] = True
                nudged = True
                loose[i]["note"] = "loose hit with at least 95 percent identity pushed strict"

        return nudged, loose



