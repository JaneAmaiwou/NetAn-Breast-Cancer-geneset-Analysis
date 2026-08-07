import os
import pandas as pd
from bioservices import KEGG

# Define the path to the output CSV file
script_dir = os.path.dirname(os.path.abspath(__file__))
output_file = os.path.abspath(os.path.join(script_dir, '../KEGG_clean.csv'))

# Print the name of the output file
print(f"Output file will be saved as: {output_file}")

# Check if the file already exists
if os.path.exists(output_file):
    print(f"The file {output_file} already exists. Exiting script.")
    exit()

# Create an instance of KEGG
k = KEGG()

# Set the organism to human
k.organism = "hsa"

# Get the list of human pathways
pathways = k.pathwayIds
print(f"Total pathways retrieved: {len(pathways)}")

# Initialize an empty list to store the data
data = []

# List of global or overview pathways to exclude
global_overview_pathways = ['01100', '01200', '01210', '01212', '01230', '01232', '01250', '01240']

# For each pathway
for i, pathway_id in enumerate(pathways):
    # Skip global or overview pathways
    if pathway_id in global_overview_pathways:
        continue

    print(f"\nProcessing pathway {i+1}/{len(pathways)}: {pathway_id}")

    # Get the pathway information
    pathway_info = k.get(pathway_id)
    print(f"Retrieved pathway information for {pathway_id}")

    # Parse the pathway information into a dictionary
    pathway_dict = k.parse(pathway_info)
    print(f"Parsed pathway information for {pathway_id}")

    # Extract the pathway name
    pathway_name = pathway_dict['NAME']

    # Check if the 'GENE' key is in the dictionary
    if 'GENE' in pathway_dict:
        # Get the list of genes for this pathway
        genes = pathway_dict['GENE']

        # For each gene
        for gene_id, gene_info in genes.items():
            # Extract the gene name
            gene_name = gene_info.split(';')[0]

            # Skip genes with no symbol or with "[KO:" or "[EC:" in their name
            if gene_name.startswith('"') or "[KO:" in gene_name or "[EC:" in gene_name:
                continue

            # Add the pathway name and gene name to the data
            data.append((pathway_name, gene_name))

print(f"\nTotal gene-pathway pairs collected: {len(data)}")

# Convert the data into a pandas DataFrame
df = pd.DataFrame(data, columns=['annotation', 'SYMBOL'])
print("\n\nDataFrame created with columns: ", df.columns)

# Convert the 'annotation' column to string type
df['annotation'] = df['annotation'].astype(str)

# Apply the replace function to clean the 'annotation' column
df['annotation'] = df['annotation'].str.replace(r"\['|- Homo sapiens \(human\)'\]", '', regex=True)

# Save the DataFrame to a CSV file
df.to_csv(output_file, index=False)
print(f"\n\nData saved to {output_file}")
