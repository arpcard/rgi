import os, argparse
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

def main(args):
    # remove dups
    records = remove_dup_seqs(SeqIO.parse(os.path.join(args.input_fasta_file), "fasta"))
    # remove sub-sequences
    final_records = remove_sub_seqs(records)
    with open(args.output_fasta_file, 'w') as fout:
        for header in final_records:
            fout.write(">{}\n{}\n".format(header, final_records[header]))

    print("Saved {} records".format(len(final_records)))

def remove_sub_seqs(records):
    matches = []
    seqs = []
    recs = {}

    for s in records:
        recs[s.id] = str(s.seq)
        if s.seq not in seqs:
            seqs.append(str(s.seq))

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if seqs[i] in seqs[j]:
                if i != j and seqs[i] != seqs[j] and seqs[i] not in matches:
                    matches.append(seqs[i])

    # remove records in matches
    final_records = {}
    for k in recs:
        if recs[k] not in matches:
           final_records[k] = recs[k] 

    return final_records

def remove_dup_seqs(records):
    """"SeqRecord iterator to removing duplicate sequences."""
    checksums = set()
    checksums_ids = set()
    for record in records:
        checksum = seguid(record.seq)
        checksum_ids = seguid(record.id)
        if checksum in checksums:# or checksum_ids in checksums_ids:
            print("Ignoring %s" % record.id)
            continue
        checksums.add(checksum)
        checksums_ids.add(checksum_ids)
        yield record

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi remove_duplicates", description='Removes duplicates sequences from annotationed fasta file')
    parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="input fasta file")
    parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="output fasta file")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()