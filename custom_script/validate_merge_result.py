import sys
from collections import defaultdict
##验证是否全部完成合并

def read_original_file(file_path):
    clusters = defaultdict(set)
    with open(file_path, 'r') as f:
        for line in f:
            gene, species, chromosome, cluster, order = line.strip().split('\t')
            clusters[cluster].add(gene)
    return clusters

def read_processed_file(file_path):
    clusters = defaultdict(set)
    chromosomes = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene, species, chromosome, cluster, order = line.strip().split('\t')
            clusters[cluster].add(gene)
            chromosomes[cluster] = chromosome
    return clusters, chromosomes

def check_naming(cluster, chromosome):
    parts = cluster.split('_')
    return parts[-1] == chromosome

def find_overlaps(clusters):
    overlaps = []
    cluster_list = list(clusters.items())
    for i in range(len(cluster_list)):
        for j in range(i+1, len(cluster_list)):
            cluster1, genes1 = cluster_list[i]
            cluster2, genes2 = cluster_list[j]
            overlap = genes1 & genes2
            if overlap and not genes1.issubset(genes2) and not genes2.issubset(genes1):
                overlaps.append((cluster1, cluster2, len(overlap)))
    return overlaps

def main(original_file, processed_file, overlap_file):
    original_clusters = read_original_file(original_file)
    processed_clusters, chromosomes = read_processed_file(processed_file)

    # Check 1: All subsets or supersets are merged
    for cluster1, genes1 in processed_clusters.items():
        for cluster2, genes2 in processed_clusters.items():
            if cluster1 != cluster2:
                if genes1.issubset(genes2) or genes2.issubset(genes1):
                    print(f"Error: Clusters {cluster1} and {cluster2} should have been merged.")

    # Check 2: Naming convention
    # for cluster, chromosome in chromosomes.items():
    #     if not check_naming(cluster, chromosome):
    #         print(f"Error: Cluster {cluster} does not follow the naming convention.")

    # Check 3: Non-overlapping original clusters are preserved
    for original_cluster, genes in original_clusters.items():
        if not any(genes.issubset(processed_genes) for processed_genes in processed_clusters.values()):
            print(f"Error: Non-overlapping original cluster {original_cluster} is missing from processed results.")

    # Check 4: All partial overlaps are in the overlap file
    processed_overlaps = find_overlaps(processed_clusters)
    with open(overlap_file, 'r') as f:
        # reported_overlaps = [tuple(line.strip().split('\t')[:2]) for line in f]
        reported_overlaps = [tuple(line.strip().split('\t')[i] for i in [0, 2]) for line in f]
    for cluster1, cluster2, _ in processed_overlaps:
        if (cluster1, cluster2) not in reported_overlaps and (cluster2, cluster1) not in reported_overlaps:
            print(f"Error: Overlap between {cluster1} and {cluster2} is missing from overlap file.")

    print("Validation complete.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate_script.py <spe>")
        sys.exit(1)

    spe = sys.argv[1]
    input_file=spe+"_genes_result"
    processed_file=spe+"_processed_genes_results"
    overlapping_file=spe+"_overlapping_clusters"
    main(input_file, processed_file, overlapping_file)