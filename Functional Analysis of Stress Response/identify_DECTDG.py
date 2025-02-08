import argparse
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats
##根据平均表达量差异与成员中差异表达基因的百分比计算差异表达基因簇
def load_data(gene_file, number_file, counts_file, anno_file):
    # 加载基因簇信息
    print("loading file ...")
    gene_clusters = pd.read_csv(gene_file, sep='\t', header=None, names=['gene_id', 'species', 'chromosome', 'cluster_id', 'gene_order'])
    cluster_numbers = pd.read_csv(number_file, sep='\t', header=None, names=['species', 'chromosome', 'cluster_id', 'gene_count'])
    
    # 加载表达数据
    counts = pd.read_csv(counts_file, sep='\t', index_col='sample')
    anno = pd.read_csv(anno_file,sep=',',index_col='sample')
    
    return gene_clusters, cluster_numbers, counts, anno

def calculate_deg(counts, anno, log2fc_threshold, padj_threshold):
    # 使用pydeseq2计算差异表达基因
    print("caculating the DEGs based on Deseq2 ...")
    dds = DeseqDataSet(
        counts=counts,
        metadata=anno,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=10,
    )
    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()
    results = stat_res.results_df

    # 筛选差异表达基因
    deg_list = results[(abs(results['log2FoldChange']) >= log2fc_threshold) & (results['padj'] < padj_threshold)].index.tolist()
    
    return deg_list, results

def calculate_deg_clusters(gene_clusters, cluster_numbers, counts, deg_list, results, anno,log2fc_threshold, padj_threshold):
    
    # 检查 counts 和 anno 的样本是否匹配
    print("caculating the CTDGs counts ...")
    counts_samples = set(counts.index)
    anno_samples = set(anno.index)
    common_samples = counts_samples & anno_samples
    
    # print("Number of samples in both counts and anno:", len(common_samples))
    # print("Samples in counts but not in anno:", counts_samples - anno_samples)
    # print("Samples in anno but not in counts:", anno_samples - counts_samples)
    # 为 gene_clusters 添加列名
    gene_clusters.columns = ['gene_id', 'species', 'chromosome', 'cluster_id', 'gene_order']
    
    # print("Number of unique genes in gene_clusters:", gene_clusters['gene_id'].nunique())
    # print("Number of unique genes in counts:", counts.columns.nunique())
    # print("Number of genes in both gene_clusters and counts:", len(set(gene_clusters['gene_id']) & set(counts.columns)))
    
    # 创建一个从基因到簇的映射
    gene_to_cluster = gene_clusters.set_index('gene_id')['cluster_id']
    
    # 只保留在 counts 中存在的基因
    gene_to_cluster = gene_to_cluster[gene_to_cluster.index.isin(counts.columns)]
    
    # print("Shape of gene_to_cluster after filtering:", gene_to_cluster.shape)
    # print("Sample of gene_to_cluster after filtering:")
    # print(gene_to_cluster.head())
    
    if gene_to_cluster.empty:
        print("Error: No genes in common between gene_clusters and counts. Cannot proceed with the analysis.")
        return []
    
    # 检查是否所有的基因都在映射中
    missing_genes = set(counts.columns) - set(gene_to_cluster.index)
    if missing_genes:
        print(f"Warning: {len(missing_genes)} genes in counts are not in gene_clusters")
    
    # 为每个基因簇计算平均表达值
    cluster_counts = counts[gene_to_cluster.index].groupby(gene_to_cluster, axis=1).mean().round().astype(int)
    
    # print("Shape of cluster_counts after groupby:", cluster_counts.shape)
    # print("Sample of cluster_counts after groupby:")
    # print(cluster_counts.head())
    
    if cluster_counts.empty:
        print("Error: cluster_counts is empty. Cannot proceed with the analysis.")
        return []
    
    
    # print("Anno columns:", anno.columns)
    # print("Anno index:", anno.index.name)
    # print("Sample of anno:")
    # print(anno.head())
    # 确保 cluster_counts 和 anno 的索引匹配
    common_samples = list(set(cluster_counts.index) & set(anno.index))
    # print("Number of common samples:", len(common_samples))
    
    if not common_samples:
        print("Error: No common samples between cluster_counts and anno. Cannot proceed with the analysis.")
        return []
    
    cluster_counts = cluster_counts.loc[common_samples]
    anno = anno.loc[common_samples]
    
    # print("Shape of cluster_counts after filtering:", cluster_counts.shape)
    # print("Shape of anno after filtering:", anno.shape)
    print("All samples in cluster_counts:")
    print(cluster_counts.index.tolist())
    print("\nAll samples in anno:")
    print(anno.index.tolist())

    # 使用DESeq2计算基因簇的差异表达
    print("caculating the DE-CTDGs based on Deseq2 ...")
    dds_cluster = DeseqDataSet(
        counts=cluster_counts,
        metadata=anno,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=10,
    )
    dds_cluster.deseq2()
    stat_res_cluster = DeseqStats(dds_cluster)
    stat_res_cluster.summary()
    results_cluster = stat_res_cluster.results_df
    
    # 筛选候选差异表达基因簇
    candidate_clusters = results_cluster[(abs(results_cluster['log2FoldChange']) >= log2fc_threshold) & (results_cluster['padj'] < padj_threshold)].index.tolist()
    
    # 确定最终的差异表达基因簇
    deg_clusters = []
    for cluster in candidate_clusters:
        cluster_genes = gene_clusters[gene_clusters['cluster_id'] == cluster]['gene_id'].tolist()
        cluster_deg_count = sum(1 for gene in cluster_genes if gene in deg_list)
        cluster_size = cluster_numbers[cluster_numbers['cluster_id'] == cluster]['gene_count'].values[0]
        if cluster_deg_count >= int(0.6 * cluster_size):
            deg_clusters.append(cluster)
    
    return deg_clusters




def main(args):
    gene_clusters, cluster_numbers, counts, anno = load_data(args.gene_file, args.number_file, args.counts_file, args.anno_file)
    deg_list, results = calculate_deg(counts, anno, args.log2fc_threshold, args.padj_threshold)
    deg_clusters = calculate_deg_clusters(gene_clusters, cluster_numbers, counts, deg_list, results,anno,args.log2fc_threshold, args.padj_threshold)
    
    # 输出差异表达基因簇列表
    with open(args.output_file, 'w') as f:
        for cluster in deg_clusters:
            f.write(f"{cluster}\n")
    print(f"all DE-CTDGs have been saved to {args.output_file}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate differentially expressed gene clusters")
    parser.add_argument("--gene_file", required=True, help="Path to the gene cluster file")
    parser.add_argument("--number_file", required=True, help="Path to the cluster number file")
    parser.add_argument("--counts_file", required=True, help="Path to the counts matrix file")
    parser.add_argument("--anno_file", required=True, help="Path to the annotation file")
    parser.add_argument("--log2fc_threshold", type=float, default=1.0, help="Log2 fold change threshold")
    parser.add_argument("--padj_threshold", type=float, default=0.05, help="Adjusted p-value threshold")
    parser.add_argument("--output_file", required=True, help="Path to the output file")
    
    args = parser.parse_args()
    main(args)
