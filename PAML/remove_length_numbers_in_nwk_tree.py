// as CODEML needs a phylogenetic tree as input, and it does not take in newick trees that have length values,
// this script is used to remove length numbers in the .nwk tree files.
import os
import re
from pathlib import Path

def remove_branch_lengths(tree_string):
    # Split the tree string using the specified separators
    tokens = re.split(r'([(),:])', tree_string)
    # Remove the numbers and colons
    tokens = [t for t in tokens if not (t.replace(".", "").isdigit() or t == ":")]
    # Reassemble the modified tree string
    return "".join(tokens)

# Set the directory containing the tree files
input_dir = "/Users/hanli/Desktop/FYP/PAML/new_fna_faa/cleaned_seq/cleaned_seq_niche_sep/aa_aligned/tree/"

# Iterate through each file in the directory
for filename in os.listdir(input_dir):
    if filename.endswith(".nwk"):  # Ensure the file has the correct extension
        filepath = os.path.join(input_dir, filename)
        
        # Read the original tree string from the file
        with open(filepath, "r") as file:
            tree_string = file.read().strip()
        
        # Remove branch lengths and save the modified tree string to a new file
        modified_tree_string = remove_branch_lengths(tree_string)
        output_filename = os.path.splitext(filename)[0] + "_nolength.nwk"
        output_filepath = os.path.join(input_dir, output_filename)
        
        with open(output_filepath, "w") as file:
            file.write(modified_tree_string)

