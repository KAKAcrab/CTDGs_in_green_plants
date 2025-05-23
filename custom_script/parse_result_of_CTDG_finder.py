import sys
from collections import defaultdict
from itertools import combinations
##对CTDG_finder结果文件genes_result.txt中存在子集或超集的基因簇进行合并
def process_input(input_file):
    cluster_to_genes = defaultdict(set)
    cluster_to_info = {}

    with open(input_file, 'r') as f:
        for line in f:
            gene, species, chromosome, cluster_name, order = line.strip().split('\t')
            cluster_to_genes[cluster_name].add(gene)
            cluster_to_info[gene] = (species, chromosome)

    return cluster_to_genes, cluster_to_info

def merge_clusters(cluster_to_genes):
    clusters_to_merge = defaultdict(set)
    
    for (name1, genes1), (name2, genes2) in combinations(cluster_to_genes.items(), 2):
        if genes1.issubset(genes2) or genes2.issubset(genes1):
            clusters_to_merge[name1].add(name2)
            clusters_to_merge[name2].add(name1)
    
    merged_clusters = {}
    processed = set()
    
    for cluster in cluster_to_genes:
        if cluster in processed:
            continue
        if cluster in clusters_to_merge:
            related_clusters = clusters_to_merge[cluster].union({cluster})
            new_name = '_'.join(sorted(related_clusters))
            new_genes = set()
            for related in related_clusters:
                new_genes.update(cluster_to_genes[related])
                processed.add(related)
            merged_clusters[new_name] = new_genes
        else:
            merged_clusters[cluster] = cluster_to_genes[cluster]
    
    return merged_clusters

def write_processed_results(merged_clusters, cluster_to_info, output_file):
    with open(output_file, 'w') as f:
        for cluster_name, genes in merged_clusters.items():
            #original_clusters = cluster_name.split('_')
            first_gene=next(iter(genes))
            species, chromosome = cluster_to_info[first_gene]#使用第一个基因的物种和染色体信息载入到基因簇中
            for i, gene in enumerate(sorted(genes), 1):
                f.write(f"{gene}\t{species}\t{chromosome}\t{cluster_name}\t{i}\n")

def find_overlapping_clusters(merged_clusters):
    overlapping = []
    for (name1, genes1), (name2, genes2) in combinations(merged_clusters.items(), 2):
        overlap = genes1.intersection(genes2)
        if overlap and not (genes1.issubset(genes2) or genes2.issubset(genes1)):
            shorter_len = min(len(genes1), len(genes2))
            overlap_ratio = len(overlap) / shorter_len
            overlapping.append((name1, len(genes1), name2, len(genes2), len(overlap), overlap_ratio, ','.join(sorted(overlap))))
    return sorted(overlapping, key=lambda x: x[5], reverse=True)

def write_overlapping_info(overlapping, output_file):
    with open(output_file, 'w') as f:
        for info in overlapping:
            f.write('\t'.join(map(str, info)) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <spe>")
        sys.exit(1)

    spe = sys.argv[1]
    input_file=spe+"_genes_result"
    processed_file=spe+"_processed_genes_results"
    overlapping_file=spe+"_overlapping_clusters"

    cluster_to_genes, cluster_to_info = process_input(input_file)
    merged_clusters = merge_clusters(cluster_to_genes)
    write_processed_results(merged_clusters, cluster_to_info, processed_file)
    overlapping = find_overlapping_clusters(merged_clusters)
    write_overlapping_info(overlapping, overlapping_file)