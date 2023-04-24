// to simplify the process and obtain some primary results from the selection pressure analysis
// this script is used to select the first 10 seqs from each gene's homolog and paralog file

import os
from Bio import SeqIO

input_dir = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa/cleaned_seq/cleaned_seq_niche_sep"
output_dir = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa/cleaned_seq/cleaned_seq_niche_sep/select_10"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file_name in os.listdir(input_dir):
    if file_name.endswith(".fasta"):
        input_file_path = os.path.join(input_dir, file_name)
        output_file_path = os.path.join(output_dir, file_name)

        with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
            seq_records = list(SeqIO.parse(input_file, "fasta"))

            if len(seq_records) > 10:
                seq_records = seq_records[:10]

            SeqIO.write(seq_records, output_file, "fasta")

