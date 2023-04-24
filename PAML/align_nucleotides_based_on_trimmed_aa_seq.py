// this script is used to perform alignment on nucleotides seqs based on aligned aa seqs.
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def back_translate(trimmed_aa_file, untrimmed_nt_file, output_file):
    trimmed_aa_sequences = SeqIO.to_dict(SeqIO.parse(trimmed_aa_file, "fasta"))
    untrimmed_nt_sequences = SeqIO.to_dict(SeqIO.parse(untrimmed_nt_file, "fasta"))

    output_sequences = []

    for record_id, aa_record in trimmed_aa_sequences.items():
        nt_sequence = untrimmed_nt_sequences[record_id].seq
        aa_sequence = aa_record.seq
        codon_aligned_sequence = []

        nt_index = 0
        for aa in aa_sequence:
            if aa == "-":
                codon_aligned_sequence.append("---")
            else:
                codon_aligned_sequence.append(nt_sequence[nt_index:nt_index+3])
                nt_index += 3

        output_sequences.append(SeqRecord(Seq("".join(codon_aligned_sequence)), id=record_id, description=""))

    SeqIO.write(output_sequences, output_file, "fasta")

input_dir = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa/cleaned_seq/cleaned_seq_niche_sep"
trim_re_aa_aligned = os.path.join(input_dir, "trim_re_aa_aligned")
cleaned_seq_niche_sep = input_dir
output_dir = os.path.join(input_dir, "trim_codon_aligned")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

trimmed_aa_files = [f for f in os.listdir(trim_re_aa_aligned) if f.endswith(".fasta")]
untrimmed_nt_files = [f for f in os.listdir(cleaned_seq_niche_sep) if f.endswith(".fasta")]

for trimmed_aa_file in trimmed_aa_files:
    KID, comparison, comparison2, order, niche = trimmed_aa_file.split("_")[3:-2]
    print(f"Current file: {trimmed_aa_file}")

    untrimmed_nt_file = [f for f in untrimmed_nt_files if f.endswith("{}_{}_{}_{}_ffn_output_{}.fasta".format(KID, comparison, comparison2, order, niche))]
    if not untrimmed_nt_file:
        print(f"Untrimmed nt file not found for {trimmed_aa_file}")
        continue
    else:
        untrimmed_nt_file = untrimmed_nt_file[0]

    output_file = "trim_codon_aligned_{}_{}_{}_{}_ffn_output_{}.fasta".format(KID, comparison, comparison2, order, niche)
    back_translate(os.path.join(trim_re_aa_aligned, trimmed_aa_file), os.path.join(cleaned_seq_niche_sep, untrimmed_nt_file), os.path.join(output_dir, output_file))

