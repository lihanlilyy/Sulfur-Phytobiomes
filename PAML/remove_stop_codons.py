# this script is used to remove stop codons in the sequences
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def remove_stop_codons_nucleotides(record):
    sequence = str(record.seq)
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    stop_codons = ['TAA', 'TAG', 'TGA']
    filtered_codons = [codon for codon in codons if codon not in stop_codons]

    new_sequence = "".join(filtered_codons)
    new_record = record
    new_record.seq = Seq(new_sequence)
    return new_record

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    records = list(SeqIO.parse(input_file, "fasta"))
    new_records = [remove_stop_codons_nucleotides(record) for record in records]

    SeqIO.write(new_records, output_file, "fasta")
