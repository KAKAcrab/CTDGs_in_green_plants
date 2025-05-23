import sys
from collections import defaultdict
from itertools import combinations
##对最终完成合并的结果中发生部分重叠的基因簇进行扫描
def read_processed_results(input_file):
    cluster_to_genes = defaultdict(set)
    with open(input_file, 'r') as f:
        for line in f:
            gene, _, _, cluster_name, _ = line.strip().split('\t')
            cluster_to_genes[cluster_name].add(gene)
    return cluster_to_genes

def find_overlapping_clusters(cluster_to_genes):
    overlapping = []
    for (name1, genes1), (name2, genes2) in combinations(cluster_to_genes.items(), 2):
        overlap = genes1.intersection(genes2)
        if overlap and not (genes1.issubset(genes2) or genes2.issubset(genes1)):
            shorter_len = min(len(genes1), len(genes2))
            overlap_ratio = len(overlap) / shorter_len
            overlapping.append((
                name1, 
                len(genes1), 
                name2, 
                len(genes2), 
                len(overlap), 
                overlap_ratio, 
                ','.join(sorted(overlap))
            ))
    return sorted(overlapping, key=lambda x: x[5], reverse=True)

def write_overlapping_info(overlapping, output_file):
    with open(output_file, 'w') as f:
        for info in overlapping:
            f.write('\t'.join(map(str, info)) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python overlap_analysis.py <spe>")
        sys.exit(1)

    spe = sys.argv[1]
    input_file=spe+"_processed_genes_results"
    out_file = spe+"_final_overlapping_clusters_info.txt"
    cluster_to_genes = read_processed_results(input_file)
    overlapping = find_overlapping_clusters(cluster_to_genes)
    write_overlapping_info(overlapping, out_file)

    print(f"Analysis complete. Results written to {out_file}")