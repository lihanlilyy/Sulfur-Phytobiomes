import os
import pandas as pd
from Bio import SeqIO

# Load data
aggregated_tests = pd.read_excel('/Users/hanli/Desktop/FYP/PAML/aggregated_all_tests_combined.xlsx', sheet_name='remove_0_rows')
final_matrix = pd.read_excel('/Users/hanli/Desktop/FYP/PAML/final_matrix_with_metadata.xlsx', sheet_name='metadata')
genes_dir = '/Users/hanli/Desktop/FYP/PAML/genomes/kofamscan/csv'
gff_dir = '/Users/hanli/Desktop/FYP/PAML/genomes/gff'
fna_dir = '/Users/hanli/Desktop/FYP/PAML/genomes/fna'

def extract_sequences(genome, gene, scaffold_id, cds_number):
    gff_file = os.path.join(gff_dir, f'{genome}.gff')
    fna_file = os.path.join(fna_dir, f'{genome}.fna')
    
    try:
        with open(gff_file) as gff:
            cds_list = []
            for line in gff:
                if line.startswith(scaffold_id) and "CDS" in line:
                    cds_list.append(line)
            cds_number_int = int(cds_number)
            cds_list = [cds for cds in cds_list if "CDS" in cds]
            if cds_number_int > len(cds_list):
                return None
            target_cds = cds_list[cds_number_int - 1]
            start, end = [int(x) for x in target_cds.split("\t")[3:5]]

            if start and end:
                with open(fna_file) as fna:
                    for record in SeqIO.parse(fna, 'fasta'):
                        if scaffold_id in record.id:
                            return str(record.seq[start - 1:end])
            else:
                return 
    except FileNotFoundError:
        print(f"File not found: {gff_file}")
        return

for _, row in aggregated_tests.iterrows():
    if row['method'] == 'KOfamScan' and row['Order'] in ['Rhizobiales', 'Pseudomonadales'] and row['selected'] == 1:
        gene = row['gene']
        niche = row['niche']
        order = row['Order']
        niche_new = niche.replace(' & ', '_')
        print("***********")
        print(gene)
        print(niche)
        print(order)
        print("***********")
        
        # Filter genomes
        if niche == "rhizo & soil":
            spheres = ["Rhizosphere", "Soil"]
        elif niche == "phyllo & soil":
            spheres = ["Phyllosphere", "Soil"]
        else:
            continue
            
        selected_genomes = final_matrix[(final_matrix['sphere'].isin(spheres)) & (final_matrix['Order'] == order)]
        
        # Extract sequences
        sequences = []
        for genome in selected_genomes['Gene']:
            csv_file = os.path.join(genes_dir, f'{genome}.csv')
            with open(csv_file) as csvfile:
                for line in csvfile:
                    if gene in line:
                        identifier, _ = line.strip().split(',', 1)
                        scaffold_id, cds_number = identifier.rsplit('_', 1)
                        scaffold_id = scaffold_id.replace('"', '')
                        cds_number = cds_number.replace('"', '')
                        seq = extract_sequences(genome, gene, scaffold_id, cds_number)
                        if seq is not None:  # Add this condition to check if the sequence is not None
                            sequences.append((genome, seq))
                        

        # Perform multiple sequence alignment using Clustal Omega
                # Save sequences to a FASTA file
        niche_order = f"{niche.replace(' & ', '_')}_{order}"
        input_file_path = f'/Users/hanli/Desktop/FYP/PAML/aligned_sequences/{gene}_{niche_order}_fna_input.fasta'
        with open(input_file_path, 'w') as fasta_file:
            for i, (genome_name, sequence) in enumerate(sequences):
                fasta_file.write(f'>{genome_name}_{i}\n{sequence}\n')
 
