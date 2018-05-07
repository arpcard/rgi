#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""
Base classes used throughout RGI
"""
import os
from Bio import SeqIO
import json
from abc import ABCMeta, abstractmethod

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

    def results(self, blast_results, query_id, perfect, strict , loose):
        """
        Sort results to perfect, strict, loose paradigm

        Parameters
        ----------

        Args:
            blast_results (dict): dictionary containing perfect, strict anf loose hits
            query_id (str): blast record query
            perfect (dict): dictionary containing perfect hits
            strict (dict): dictionary containing strict hits
            loose (dict): dictionary containing loose hits

        Returns:
            blast_results (dict): dictionary of sorted results
        """     
        if len(perfect) == 0 and len(strict) == 0 and len(loose) > 0:
            nudged , loose = self.nudge(loose)
            if nudged is True and self.loose is False:
                blast_results[query_id] = loose
            elif self.loose is True:
                blast_results[query_id] = loose

        elif len(perfect) == 0:
            if len(strict) > 0: 
                blast_results[query_id] = strict
        else:
            if len(perfect) > 0: 
                blast_results[query_id] = perfect

        return blast_results

    def nudge(self, loose):
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
            # half_pass_bitscore = int(float(loose[i]["pass_bitscore"])/2)
            # if 95 <= int(loose[i]["perc_identity"]) <= 100 and loose[i]["bit_score"] > half_pass_bitscore:
            if 95 <= int(loose[i]["perc_identity"]) <= 100:
                # add to strict 
                loose[i]["type_match"] = "Strict"
                nudged = True
                # print("Nudged Loose -> Strict ... HSP: {} | model_name: {} | perc_identity: {}% |  bit_score : {} | pass_bitscore:  {} | half_pass_bitscore : {}"\
                #     .format(i,loose[i]["model_name"], loose[i]["perc_identity"], loose[i]["bit_score"], loose[i]["pass_bitscore"], half_pass_bitscore))
        return nudged, loose



