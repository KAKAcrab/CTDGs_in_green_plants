import sys
import pandas as pd

def count_gene_clusters(input_file, output_file):
    # Read the input file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['GeneID', 'SpeciesID', 'ChromosomeID', 'ClusterID', 'Order'])

    # Group by SpeciesID, ChromosomeID, and ClusterID, and count the number of genes in each cluster
    cluster_counts = df.groupby(['SpeciesID', 'ChromosomeID', 'ClusterID']).size().reset_index(name='MemberCount')

    # Save the result to an output file
    cluster_counts.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        count_gene_clusters(input_file, output_file)
