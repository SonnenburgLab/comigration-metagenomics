import os
import sys
import glob
import pandas as pd
from Bio import SeqIO

input_cluster_file = sys.argv[1]
genome_dir = sys.argv[2]
species = sys.argv[3]
ref_genome = sys.argv[4]
ref_sample = ref_genome.split('__')[1]

print(sys.argv)


def count_files_with_extension(directory, extension):
    # Construct the pattern to match files with the specified extension
    pattern = os.path.join(directory, f'*.{extension}')
    # Use glob to find files matching the pattern
    matching_files = glob.glob(pattern) 
    return len(matching_files)


def read_contig_names(fasta_file):
    """
    Reads contig names from a FASTA file and returns them as a list.

    Args:
    - fasta_file (str): Path to the FASTA file.

    Returns:
    - list: List of contig names.
    """
    contig_names = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig_names.append(record.id)
    return contig_names

ref_genome_gene_names = read_contig_names(os.path.join('species', species, 'prodigal', ref_genome + '.faa'))

df = pd.read_csv(input_cluster_file, sep='\t', header=None)

df.rename(columns={0: "Target", 1: "Query"}, inplace=True)

grouped_df = df.groupby('Target').size()

threshold = count_files_with_extension(genome_dir, 'fa')*0.9

core_genes = grouped_df[grouped_df > threshold].index.tolist()

core_filtered_clusters = df[df['Target'].isin(core_genes)]
print(core_filtered_clusters)

# Define the string you want to match
match_string = '%s__' % ref_sample  # Change this to your desired string

# Filter the dataframe
ref_genome_filtered_df = core_filtered_clusters[core_filtered_clusters.apply(lambda row: any(match_string in str(cell.replace('-', '_')) for cell in row), axis=1)]

ref_genome_core_genes = ref_genome_filtered_df['Query'].tolist()

ref_genome_core_genes = list(set(core_filtered_clusters['Query'].tolist()).intersection(set(ref_genome_gene_names)))

# Define the filename
out_filename = os.path.join("species", species, "mmseqs", species+"_core_cluster.tsv")

print(out_filename)

# Write the categories to the file
with open(out_filename, 'w') as file:
    for gene in ref_genome_core_genes:
        file.write(gene + '\n')
