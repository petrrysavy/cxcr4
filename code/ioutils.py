from Bio import SeqIO


def load_single_fasta_sequence(path):
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq)
    raise RuntimeError("CXCR4 sequence file does not contain the sequence")