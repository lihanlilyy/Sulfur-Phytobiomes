# this script is used to remove the extra seqs that do not correspond to the target gene.
# the target gene names are given in the keywords list below.
import os
import re
from Bio import SeqIO

input_directory = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa"
output_directory = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa/cleaned_seq"

keywords = {
    "K00548": "Methionine synthase",
    "K00651": "O-acetyltransferase",
    "K00799": "S-transferase",
    "K01505": "deaminase",
    "K07160": "5-oxoprolinase",
    "K15554": "permease",
    "K22955": "lactone",
    "K00003": "dehydrogenase",
    "K00032": "reductase",
    "K00036": "dehydrogenase",
    "K00383": "reductase",
    "K02439": "sulfurtransferase",
    "K00387": "sulfoxide",
    "K00432": "gpx|btuE|bsaA",
    "K00548": "synthase",
    "K00640": "acetyltransferase",
    "K05301": "reductase",
    "K06048": "ligase",
    "K15554": "permease",
    "K17228": "monooxygenase",
}

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for file_name in os.listdir(input_directory):
    if file_name.startswith(tuple(keywords.keys())) and file_name.endswith(".fasta"):
        input_file_path = os.path.join(input_directory, file_name)
        output_file_path = os.path.join(output_directory, file_name)

        with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file:
            fasta_records = list(SeqIO.parse(input_file, "fasta"))

            for record in fasta_records:
                key = file_name.split("_")[0]

                if re.search(keywords[key], record.description, re.IGNORECASE):
                    SeqIO.write(record, output_file, "fasta")

