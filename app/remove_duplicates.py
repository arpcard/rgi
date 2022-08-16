import os, argparse
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from app.settings import *
from collections import defaultdict
from argparse import RawTextHelpFormatter
from app.settings import APP_NAME, SOFTWARE_VERSION

def main(args):
    if args.debug:
        logger.setLevel(10)

    dedup_records = defaultdict(list)
    card_cannonical_lengths = {}
    # find the lengths of all reference cannonical sequences
    for record in SeqIO.parse(os.path.join(args.card_annotation), "fasta"):
        if record.id[:4] == "ARO:":
            accession = record.id.split("|")[0].split(":")[1]
            card_cannonical_lengths[accession] = len(record.seq)

    remove = []
    # use reference cannonical sequences lengths to filter out shorter variants
    for record in SeqIO.parse(os.path.join(args.input_fasta_file), "fasta"):
        if record.id[:23] == "Prevalence_Sequence_ID:":
            accession = record.id.split("|")[3].split(":")[1].strip("")
            cl = card_cannonical_lengths[accession]
            pl = len(record.seq)
            perc = int((pl/cl)*100)
            # if cl > pl and (int(perc) < 90):
            if cl != pl:
                logger.info("remove: {} -- {} > {} ({})".format(record.id, cl, pl, perc))
                if record.id not in remove:
                    remove.append(record.id)
            else:
                dedup_records[str(record.seq)].append(record.id)

    # write FASTA file
    logger.info("write FASTA file ...")
    with open(args.output_fasta_file, 'w') as fout:
        for seq in dedup_records:
            fout.write(">{}\n{}\n".format(dedup_records[seq][-1], seq))

    logger.info("Number of card cannonical records: {}".format(len(card_cannonical_lengths.keys())))
    logger.info("Number of prevalence sequences removed based on length of reference: {}".format(len(remove)))
    logger.info("Saved {} records".format(len(dedup_records)))
    logger.info("Done.")

def remove_sub_sequences(records):
    """
    Remove sub sequences. i.e sequence contained within other sequences

    Parameters
    ----------

    Args:
        records (dict): dictionary containing sequences

    Returns:
        records (dict): dictionary containing sequences
    """
    matches = []
    seqs = []
    recs = {}

    for s in records:
        recs[s.id] = str(s.seq)
        if s.seq not in seqs:
            seqs.append(str(s.seq))

    logger.info("1: {}".format(len(recs.keys())))
    logger.info("2: {}".format(len(seqs)))

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if seqs[i] in seqs[j]:
                if i != j and seqs[i] != seqs[j] and seqs[i] not in matches:
                    matches.append(seqs[i])

    logger.info("3: {}".format(len(matches)))
    # remove records in matches
    final_records = {}
    for k in recs:
        if recs[k] not in matches:
           final_records[k] = recs[k]

    return final_records


def remove_duplicate_sequences(records):
    """"SeqRecord iterator to removing duplicate sequences."""
    checksums = set()
    all_seqs = []
    for record in records:
        all_seqs.append(record.seq)
        checksum = seguid(record.seq)
        if checksum in checksums:
            logger.warning("Ignoring {}".format(record.id))
            continue
        checksums.add(checksum)

        yield record

    logger.info("Number of sequences: {}".format(len(all_seqs)))

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi remove_duplicates", description='{} - {} - Remove duplicates \n\nRemoves duplicates sequences from annotationed fasta file.'.format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="input fasta file")
    parser.add_argument('--card_annotation', dest="card_annotation", required=True, help="card_annotation input fasta file")
    parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="output fasta file")
    parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
