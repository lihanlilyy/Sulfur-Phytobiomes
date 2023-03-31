import pandas as pd
from ete3 import Tree

def extract_genome_name(node_name):
    parts = node_name.split("_")
    if parts[0].startswith("Ga"):
        return parts[0]
    else:
        return parts[0] +"_" +parts[1]

# Read metadata from Excel file
metadata_file = "/Users/hanli/Desktop/FYP/sulfer_genes/final_matrix_with_metadata.xlsx"
metadata = pd.read_excel(metadata_file, sheet_name="metadata")

# Create a dictionary to map genome names to their sphere values
genome_to_sphere = {}
for index, row in metadata.iterrows():
    genome_to_sphere[row["Gene"]] = row["sphere"]

# Load the tree file
tree_file = "/Users/hanli/Downloads/clustal-omega-1.2.4/src/outputs/trees/K00432_phyllo_soil_Rhizobiales_fna_input.tree"
tree = Tree(tree_file)

# Define the colors for the different sphere values
sphere_colors = {
    "Rhizosphere": "#FFD700",
    "Phyllosphere": "#3CB371",
    "Soil": "#B8860B",
}

# Generate iTOL dataset file for peripheral ring colors
# itol_ring_colors_file = "itol_ring_colors.txt"
# with open(itol_ring_colors_file, "w") as itol_file:
#     itol_file.write("DATASET_COLORSTRIP\n")
#     itol_file.write("SEPARATOR TAB\n")
#     itol_file.write("DATASET_LABEL\tsphere_colors\n")
#     itol_file.write("COLOR\t#ff0000\n")
#     itol_file.write("DATA\n")
    
#     for leaf in tree.iter_leaves():
#         genome_name = extract_genome_name(leaf.name)
#         print(genome_name)
#         sphere = genome_to_sphere.get(genome_name, None)
#         if sphere:
#             color = sphere_colors.get(sphere, "#000000")  # Default to black if the sphere is not found
#             itol_file.write(f"{leaf.name}\t{color}\n")

# Generate iTOL dataset file for branch colors
itol_branch_colors_file = "itol_branch_colors.txt"
with open(itol_branch_colors_file, "w") as itol_file:
    itol_file.write("TREE_COLORS\n")
    itol_file.write("SEPARATOR TAB\n")
    itol_file.write("DATASET_LABEL\tniche_colors\n")
    itol_file.write("COLOR\t#ff0000\n")
    itol_file.write("DATA\n")

    for leaf in tree.iter_leaves():
        color = ""
        if "Rhizobiales" in leaf.name:
            color = "#FF0000"
        elif "Pseudomonadales" in leaf.name:
            color = "#0000FF"
        
        if color:
            itol_file.write(f"{leaf.name}\tlabel\t{color}\n")

