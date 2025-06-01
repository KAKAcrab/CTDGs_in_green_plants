#!/usr/bin/env python3
"""
DE-CTDG-v2: 差异表达串联重复基因簇(DE-CTDGs)分析工具
作者: yyh
日期: 2025-04-12

基于原始DE_CTDG_for_pear.py脚本，新增分析成员数量为2的DE-CTDGs中两个成员的表达差异
   - 提取所有两个成员的DE-CTDGs在各种条件下的表达量
   - 根据表达量高低，将成员分为"major"和"minor"
   - 输出表格，包含基因簇ID以及"major"和"minor"成员的表达量

满足以下条件即认定为DE-CTDGs：
1. 任一基因簇80%的成员（省去计算得到的小数部分，只保留整数）均为差异表达基因（DEGs，基于deseq2算法）
2. 使用所有成员的平均表达量作为该CTDG的表达量并将CTDG视为一个基因，然后满足deseq2算法的差异表达基因

用法:
    python enhanced_dectdg_analyzer.py \
      --input_dir /path/to/data \
      --gff_file /path/to/annotation.gff \
      --gene_file /path/to/processed_raw_genes_result \
      --number_file /path/to/processed_raw_numbers_result \
      --output_dir /path/to/output \
      --threads 8
"""

import pandas as pd
import numpy as np
import argparse
import os
import logging
import time
import sys
import re
import warnings
import random
from itertools import groupby
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
from scipy import stats
from collections import defaultdict
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 设置日志格式
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', 
                    level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_gff(gff_file):
    """
    从GFF文件中提取mRNA的ID和长度信息
    只提取第三列为"mRNA"且第九列"ID="对应的基因ID信息
    
    参数:
        gff_file: GFF文件路径
        
    返回:
        gene_lengths: 基因ID和长度的字典
    """
    logger.info(f"解析GFF文件: {gff_file}")
    gene_lengths = {}
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                # 只处理第三列为"mRNA"的行
                if parts[2] == "gene":
                # if parts[2] == "mRNA":
                    # 提取属性
                    attributes = parts[8]
                    gene_id = None
                    
                    # 仅使用ID=模式提取基因ID
                    match = re.search(r'ID=([^;]+)', attributes)
                    if match:
                        gene_id = match.group(1)
                    
                    if gene_id:
                        # 计算基因长度
                        try:
                            start = int(parts[3])
                            end = int(parts[4])
                            length = end - start + 1
                            gene_lengths[gene_id] = length
                        except ValueError:
                            logger.warning(f"无法解析坐标: {parts[3]}-{parts[4]}")
    except Exception as e:
        logger.error(f"解析GFF文件时出错: {str(e)}")
        raise
    
    logger.info(f"从GFF文件中提取了 {len(gene_lengths)} 个mRNA的ID和长度信息")
    return gene_lengths

def load_gene_clusters(gene_file, number_file):
    """
    加载基因簇信息
    
    参数:
        gene_file: 基因簇文件路径 (_processed_raw_genes_result)
        number_file: 基因簇大小文件路径 (_processed_raw_numbers_result)
        
    返回:
        gene_clusters: 基因到簇的映射DataFrame
        cluster_numbers: 簇大小信息DataFrame
    """
    logger.info(f"加载基因簇信息: {gene_file}")
    try:
        # 检查第一行是否为标题行
        with open(gene_file, 'r') as f:
            first_line = f.readline().strip().split('\t')
            if first_line and len(first_line) >= 5 and 'gene' in first_line[0].lower():
                logger.info(f"检测到标题行，跳过第一行: {first_line}")
                skip_rows = 1
            else:
                skip_rows = 0
        
        gene_clusters = pd.read_csv(gene_file, sep='\t', header=None, skiprows=skip_rows,
                                  names=['gene_id', 'species', 'chromosome', 'cluster_id', 'gene_order'])
    except Exception as e:
        logger.error(f"加载基因簇文件时出错: {str(e)}")
        raise
    
    logger.info(f"加载基因簇大小信息: {number_file}")
    try:
        # 检查第一行是否为标题行
        with open(number_file, 'r') as f:
            first_line = f.readline().strip().split('\t')
            if first_line and len(first_line) >= 4 and 'member' in first_line[3].lower():
                logger.info(f"检测到标题行，跳过第一行: {first_line}")
                skip_rows = 1
            else:
                skip_rows = 0
        
        cluster_numbers = pd.read_csv(number_file, sep='\t', header=None, skiprows=skip_rows,
                                    names=['species', 'chromosome', 'cluster_id', 'gene_count'])
        
        # 确保gene_count列为整数
        try:
            cluster_numbers['gene_count'] = cluster_numbers['gene_count'].astype(int)
            logger.info(f"成功将gene_count列转换为整数类型")
        except (ValueError, TypeError) as e:
            logger.warning(f"无法将gene_count列转换为整数，可能包含非数字值: {str(e)}")
            # 尝试逐行转换，保留可以转换的值
            valid_rows = []
            for idx, row in cluster_numbers.iterrows():
                try:
                    gene_count = int(row['gene_count'])
                    # 更新值
                    row['gene_count'] = gene_count
                    valid_rows.append(row)
                except (ValueError, TypeError):
                    logger.warning(f"跳过无效的gene_count行: {row['cluster_id']} - {row['gene_count']}")
            
            # 重建DataFrame只包含有效行
            if valid_rows:
                cluster_numbers = pd.DataFrame(valid_rows)
                logger.info(f"清理后保留了 {len(valid_rows)} 行有效数据")
            else:
                logger.error("没有找到有效的gene_count数据")
                raise ValueError("No valid gene_count data found")
    except Exception as e:
        logger.error(f"加载基因簇大小文件时出错: {str(e)}")
        raise
    
    # 检查数据完整性
    if gene_clusters.empty or cluster_numbers.empty:
        logger.error("基因簇数据为空")
        raise ValueError("Empty gene cluster data")
    
    logger.info(f"加载了 {len(gene_clusters)} 条基因簇成员数据和 {len(cluster_numbers)} 条基因簇大小数据")
    return gene_clusters, cluster_numbers

def load_expression_data(counts_file, anno_file):
    """
    加载表达量数据和样本注释信息
    
    参数:
        counts_file: 表达量矩阵文件路径 (_counts_matrix.tsv)
        anno_file: 样本注释文件路径 (_pasAnno.csv)
        
    返回:
        counts: 表达量矩阵DataFrame
        anno: 样本注释DataFrame
    """
    logger.info(f"加载表达量矩阵: {counts_file}")
    try:
        # 检查第一行是否为标题行
        with open(counts_file, 'r') as f:
            first_line = f.readline().strip()
            if first_line and not first_line.startswith('sample') and 'GeneID' in first_line:
                # 文件有标题但没有明确的索引列
                counts = pd.read_csv(counts_file, sep='\t')
                # 如果第一列没有列名，设置为'sample'
                if counts.columns[0] == '0' or counts.columns[0] == 'Unnamed: 0':
                    counts = pd.read_csv(counts_file, sep='\t', index_col=0)
                    counts.index.name = 'sample'
                else:
                    # 将第一列设置为索引
                    first_col = counts.columns[0]
                    counts = counts.set_index(first_col)
            else:
                # 标准格式，带索引列名
                counts = pd.read_csv(counts_file, sep='\t', index_col=0)
    except Exception as e:
        logger.error(f"加载表达量矩阵时出错: {str(e)}")
        if logger.level <= logging.DEBUG:
            import traceback
            logger.debug(f"错误详细信息: {traceback.format_exc()}")
        raise
    
    logger.info(f"加载样本注释: {anno_file}")
    # 先尝试加载原始路径
    try:
        # 先检查文件是否存在
        if not os.path.exists(anno_file):
            # 尝试tsv后缀
            if anno_file.endswith('.csv'):
                alt_anno_file = anno_file.replace('.csv', '.tsv')
                if os.path.exists(alt_anno_file):
                    logger.info(f"找到替代注释文件: {alt_anno_file}")
                    anno_file = alt_anno_file
                else:
                    logger.error(f"未找到注释文件: {anno_file}")
                    raise FileNotFoundError(f"未找到注释文件: {anno_file}")
        
        # 先检查分隔符
        with open(anno_file, 'r') as f:
            first_line = f.readline().strip()
            if ',' in first_line:
                sep = ','
            elif '\t' in first_line:
                sep = '\t'
            else:
                sep = ','  # 默认为CSV
                
        # 检查是否有标题行
        with open(anno_file, 'r') as f:
            first_line = f.readline().strip().split(sep)
            if first_line and len(first_line) >= 2 and 'sample' in first_line[0].lower():
                # 有标题行，指定索引列
                anno = pd.read_csv(anno_file, sep=sep, index_col=0)
            else:
                # 无标题行或标题行不标准
                try:
                    # 尝试按第一列作为索引读取
                    anno = pd.read_csv(anno_file, sep=sep, header=None, 
                                     names=['sample', 'condition'], index_col=0)
                except Exception:
                    # 如果失败，不指定索引列读取然后手动设置索引
                    anno = pd.read_csv(anno_file, sep=sep, header=None, 
                                     names=['sample', 'condition'])
                    anno = anno.set_index('sample')
    except Exception as e:
        logger.error(f"加载样本注释时出错: {str(e)}")
        if logger.level <= logging.DEBUG:
            import traceback
            logger.debug(f"错误详细信息: {traceback.format_exc()}")
        raise
    
    # 确保样本名称匹配
    common_samples = list(set(counts.index) & set(anno.index))
    if not common_samples:
        logger.error("表达量矩阵和样本注释中没有匹配的样本")
        logger.debug(f"counts样本: {counts.index.tolist()}")
        logger.debug(f"anno样本: {anno.index.tolist()}")
        raise ValueError("No common samples between counts and annotation")
    
    counts = counts.loc[common_samples]
    anno = anno.loc[common_samples]
    
    # 检查注释列
    if 'condition' not in anno.columns:
        logger.error("样本注释中缺少'condition'列")
        logger.debug(f"anno列: {anno.columns.tolist()}")
        raise ValueError("Missing 'condition' column in annotation file")
    
    # 检查处理条件
    conditions = anno['condition'].unique()
    if 'treated' not in conditions or 'untreated' not in conditions:
        logger.warning(f"样本注释中的条件不是标准的'treated'和'untreated'，实际值为: {conditions}")
    
    logger.info(f"共找到 {len(common_samples)} 个匹配的样本")
    return counts, anno

def calculate_tpm(counts, gene_lengths):
    """
    计算TPM (Transcripts Per Million) 值
    
    参数:
        counts: 表达量矩阵DataFrame
        gene_lengths: 基因长度字典
        
    返回:
        tpm: TPM值矩阵DataFrame
    """
    logger.info("计算TPM值")
    
    # 创建新的DataFrame来存储TPM值
    tpm = pd.DataFrame(index=counts.index, columns=counts.columns)
    
    # 只计算在gene_lengths中有长度信息的基因
    valid_genes = [gene for gene in counts.columns if gene in gene_lengths]
    
    if not valid_genes:
        logger.error("没有找到任何有长度信息的基因")
        raise ValueError("No genes with length information")
    
    logger.info(f"计算 {len(valid_genes)} 个基因的TPM值")
    
    # 对每个样本计算TPM
    for sample in counts.index:
        # 首先计算RPK (Reads Per Kilobase)
        rpk_values = {}
        for gene in valid_genes:
            try:
                rpk_values[gene] = counts.loc[sample, gene] / (gene_lengths[gene] / 1000)
            except ZeroDivisionError:
                logger.warning(f"基因 {gene} 的长度为0，设置为最小值1")
                rpk_values[gene] = counts.loc[sample, gene] / 0.001  # 假设最小长度为1bp
            except KeyError:
                continue
                
        # 计算每百万的乘数
        per_million = sum(rpk_values.values()) / 1e6
        
        if per_million > 0:
            # 计算TPM
            for gene in rpk_values:
                tpm.loc[sample, gene] = rpk_values[gene] / per_million
        else:
            logger.warning(f"样本 {sample} 的总RPK为0，TPM值将全部为0")
    
    # 对于没有计算的基因，填充为0
    tpm = tpm.fillna(0)
    
    return tpm

def process_condition(condition_data):
    """
    处理单个条件的差异表达分析
    
    参数:
        condition_data: 包含条件所需数据的字典
        
    返回:
        包含条件处理结果的字典
    """
    counts = condition_data['counts']
    anno = condition_data['anno']
    gene_clusters = condition_data['gene_clusters']
    cluster_numbers = condition_data['cluster_numbers']
    log2fc_threshold = condition_data['log2fc_threshold']
    padj_threshold = condition_data['padj_threshold']
    tpm = condition_data['tpm']
    condition_name = condition_data['condition_name']
    
    logger.info(f"处理条件: {condition_name}")
    
    try:
        # 1. 使用PyDESeq2计算差异表达基因
        logger.info(f"计算 {condition_name} 条件下的差异表达基因")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            dds = DeseqDataSet(
                counts=counts,
                metadata=anno,
                design_factors="condition",
                refit_cooks=True,
                n_cpus=condition_data['threads_per_condition'],
            )
            dds.deseq2()
            stat_res = DeseqStats(dds)
            stat_res.summary()
            results = stat_res.results_df
        
        # 筛选差异表达基因
        deg_list = results[(abs(results['log2FoldChange']) >= log2fc_threshold) & 
                          (results['padj'] < padj_threshold)].index.tolist()
        logger.info(f"在 {condition_name} 条件下找到 {len(deg_list)} 个差异表达基因")
        
        # 2. 计算基因簇的平均表达值
        gene_to_cluster = gene_clusters.set_index('gene_id')['cluster_id']
        gene_to_cluster = gene_to_cluster[gene_to_cluster.index.isin(counts.columns)]
        
        if gene_to_cluster.empty:
            logger.error(f"{condition_name} 条件下没有找到匹配的基因")
            raise ValueError(f"No matching genes for condition {condition_name}")
        
        logger.info(f"计算 {condition_name} 条件下基因簇的平均表达值")
        cluster_counts = counts[gene_to_cluster.index].groupby(gene_to_cluster, axis=1).mean().round().astype(int)
        
        # 3. 计算基因簇的TPM值
        logger.info(f"计算 {condition_name} 条件下基因簇的TPM值")
        cluster_tpm = tpm[gene_to_cluster.index].groupby(gene_to_cluster, axis=1).mean()
        
        # 4. 计算基因簇的差异表达
        logger.info(f"计算 {condition_name} 条件下的差异表达基因簇")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            dds_cluster = DeseqDataSet(
                counts=cluster_counts,
                metadata=anno,
                design_factors="condition",
                refit_cooks=True,
                n_cpus=condition_data['threads_per_condition'],
            )
            dds_cluster.deseq2()
            stat_res_cluster = DeseqStats(dds_cluster)
            stat_res_cluster.summary()
            results_cluster = stat_res_cluster.results_df
        
        # 筛选候选差异表达基因簇
        candidate_clusters = results_cluster[(abs(results_cluster['log2FoldChange']) >= log2fc_threshold) & 
                                          (results_cluster['padj'] < padj_threshold)].index.tolist()
        logger.info(f"在 {condition_name} 条件下找到 {len(candidate_clusters)} 个候选差异表达基因簇")
        
        # 5. 确定最终的差异表达基因簇
        deg_clusters = []
        for cluster in candidate_clusters:
            try:
                cluster_genes = gene_clusters[gene_clusters['cluster_id'] == cluster]['gene_id'].tolist()
                cluster_deg_count = sum(1 for gene in cluster_genes if gene in deg_list)
                
                cluster_size_rows = cluster_numbers[cluster_numbers['cluster_id'] == cluster]
                if not cluster_size_rows.empty:
                    # 确保 cluster_size 是整数
                    try:
                        cluster_size = int(cluster_size_rows.iloc[0]['gene_count'])
                        # 使用80%的阈值，取整数部分
                        threshold = int(0.8 * cluster_size)
                        if cluster_deg_count >= threshold:
                            deg_clusters.append(cluster)
                    except (ValueError, TypeError):
                        logger.warning(f"基因簇 {cluster} 的大小格式不正确: {cluster_size_rows.iloc[0]['gene_count']}")
            except Exception as e:
                logger.warning(f"处理基因簇 {cluster} 时出错: {str(e)}")
        
        logger.info(f"在 {condition_name} 条件下找到 {len(deg_clusters)} 个差异表达基因簇")
        
        # 获取处理组样本
        treated_samples = anno[anno['condition'] == 'treated'].index.tolist()
        
        # 计算处理组的平均TPM
        treated_tpm = cluster_tpm.loc[treated_samples].mean() if treated_samples else pd.Series()
        
        # 存储每个基因簇的成员基因表达量 (用于后续分析成员数量为2的DE-CTDGs)
        member_expression = {}
        
        # 仅处理成员数量为2的基因簇
        for cluster in deg_clusters:
            # 获取该簇的基因数量
            cluster_size_rows = cluster_numbers[cluster_numbers['cluster_id'] == cluster]
            if not cluster_size_rows.empty:
                try:
                    cluster_size = int(cluster_size_rows.iloc[0]['gene_count'])
                    if cluster_size == 2:
                        # 获取该簇的所有成员基因
                        cluster_genes = gene_clusters[gene_clusters['cluster_id'] == cluster]['gene_id'].tolist()
                        
                        # 仅保留同时存在于TPM矩阵中的基因
                        valid_genes = [gene for gene in cluster_genes if gene in tpm.columns]
                        
                        if len(valid_genes) == 2:  # 确保有两个有效基因
                            # 计算处理组的平均TPM
                            gene1_tpm = tpm.loc[treated_samples, valid_genes[0]].mean()
                            gene2_tpm = tpm.loc[treated_samples, valid_genes[1]].mean()
                            
                            # 确定major和minor
                            if gene1_tpm >= gene2_tpm:
                                major_gene = valid_genes[0]
                                minor_gene = valid_genes[1]
                                major_tpm = gene1_tpm
                                minor_tpm = gene2_tpm
                            else:
                                major_gene = valid_genes[1]
                                minor_gene = valid_genes[0]
                                major_tpm = gene2_tpm
                                minor_tpm = gene1_tpm
                            
                            # 存储结果
                            member_expression[cluster] = {
                                'major_gene': major_gene,
                                'minor_gene': minor_gene,
                                'major_tpm': major_tpm,
                                'minor_tpm': minor_tpm
                            }
                except (ValueError, TypeError):
                    logger.warning(f"处理成员数量为2的基因簇 {cluster} 时出错")
        
        # 返回处理结果
        return {
            'condition_name': condition_name,
            'deg_list': deg_list,
            'deg_clusters': deg_clusters,
            'cluster_tpm': cluster_tpm,
            'treated_tpm': treated_tpm,  # 处理组的平均TPM
            'member_expression': member_expression,  # 成员数量为2的DE-CTDGs的成员表达量
            'results_cluster': results_cluster  # 保存基因簇的差异表达结果
        }
    
    except Exception as e:
        logger.error(f"处理条件 {condition_name} 时出错: {str(e)}")
        raise

def analyze_member_expression(condition_results):
    """
    分析成员数量为2的DE-CTDGs中两个成员的表达差异
    
    参数:
        condition_results: 条件处理结果列表
        
    返回:
        member_expression_results: 按条件分组的成员表达分析结果
    """
    logger.info("分析成员数量为2的DE-CTDGs中两个成员的表达差异")
    
    member_expression_results = {}
    
    for result in condition_results:
        condition_name = result['condition_name']
        member_expression = result['member_expression']
        
        if not member_expression:
            logger.warning(f"条件 {condition_name} 下没有成员数量为2的DE-CTDGs")
            continue
        
        # 创建该条件下的结果DataFrame
        member_df = pd.DataFrame(columns=['cluster_id', 'major_gene', 'minor_gene', 'major_tpm', 'minor_tpm', 'ratio'])
        
        # 填充数据
        for cluster, data in member_expression.items():
            # 计算major与minor的比例
            ratio = data['major_tpm'] / data['minor_tpm'] if data['minor_tpm'] > 0 else float('inf')
            
            # 添加到DataFrame
            member_df = pd.concat([member_df, pd.DataFrame({
                'cluster_id': [cluster],
                'major_gene': [data['major_gene']],
                'minor_gene': [data['minor_gene']],
                'major_tpm': [data['major_tpm']],
                'minor_tpm': [data['minor_tpm']],
                'ratio': [ratio]
            })], ignore_index=True)
        
        # 按major_tpm降序排序
        member_df = member_df.sort_values(by='major_tpm', ascending=False)
        
        # 存储结果
        member_expression_results[condition_name] = member_df
        
        logger.info(f"完成条件 {condition_name} 下 {len(member_df)} 个成员数量为2的DE-CTDGs的表达分析")
    
    return member_expression_results

def create_expression_matrices(all_clusters, condition_results, cluster_numbers):
    """
    创建表达矩阵
    
    参数:
        all_clusters: 所有基因簇ID列表
        condition_results: 不同条件的处理结果列表
        cluster_numbers: 簇大小信息DataFrame
        
    返回:
        binary_matrix: 01分布矩阵
        tpm_matrix: TPM表达量矩阵
    """
    logger.info("创建表达矩阵")
    
    try:
        # 1. 创建01矩阵
        conditions = [result['condition_name'] for result in condition_results]
        binary_matrix = pd.DataFrame(0, index=all_clusters, columns=conditions)
        
        # 添加size列
        size_mapping = {}
        
        # 只处理cluster_id存在于all_clusters中的行
        valid_clusters = cluster_numbers[cluster_numbers['cluster_id'].isin(all_clusters)]
        logger.info(f"用于创建表达矩阵的有效基因簇: {len(valid_clusters)}")
        
        for _, row in valid_clusters.iterrows():
            try:
                cluster_id = row['cluster_id']
                # 确保size是整数
                size = int(row['gene_count'])
                
                if size == 2:
                    size_mapping[cluster_id] = "member=2"
                elif 2 < size <= 5:
                    size_mapping[cluster_id] = "3<member<5"
                elif 5 < size <= 10:
                    size_mapping[cluster_id] = "6<member<10"
                else:
                    size_mapping[cluster_id] = "member>10"
            except (ValueError, TypeError) as e:
                logger.warning(f"转换基因簇 {row['cluster_id']} 的大小时出错: {str(e)}")
                size_mapping[row['cluster_id']] = "unknown"
        
        # 对没有在size_mapping中的基因簇设置为unknown
        binary_matrix['size'] = binary_matrix.index.map(lambda x: size_mapping.get(x, "unknown"))
        
        # 填充01矩阵
        for result in condition_results:
            condition_name = result['condition_name']
            deg_clusters = result['deg_clusters']
            
            if deg_clusters:  # 确保列表不为空
                binary_matrix.loc[binary_matrix.index.isin(deg_clusters), condition_name] = 1
        
        # 重新排列列，使得'size'列在最前面
        cols = binary_matrix.columns.tolist()
        size_cols = ['size']
        data_cols = [col for col in cols if col != 'size']
        binary_matrix = binary_matrix[size_cols + data_cols]
        
        # 2. 创建TPM矩阵
        tpm_matrix = pd.DataFrame(index=all_clusters)
        
        # 对每个条件，添加log2(TPM+1)值
        for result in condition_results:
            condition_name = result['condition_name']
            treated_tpm = result['treated_tpm']
            
            # 应用log2(TPM+1)转换
            tpm_log2 = np.log2(treated_tpm + 1)
            
            # 将结果添加到矩阵中
            for cluster in all_clusters:
                if cluster in tpm_log2:
                    tpm_matrix.loc[cluster, condition_name] = tpm_log2[cluster]
                else:
                    tpm_matrix.loc[cluster, condition_name] = 0
        
        return binary_matrix, tpm_matrix
    except Exception as e:
        logger.error(f"创建表达矩阵时出错: {str(e)}")
        if logger.level <= logging.DEBUG:
            import traceback
            logger.debug(f"错误详细信息: {traceback.format_exc()}")
        raise

def get_all_conditions(input_dir):
    """
    从输入目录获取所有处理条件
    
    参数:
        input_dir: 输入文件目录
        
    返回:
        处理条件列表
    """
    conditions = []
    pattern = r'(.+)_counts_matrix\.tsv'
    
    try:
        files = os.listdir(input_dir)
        logger.debug(f"目录 {input_dir} 中的文件: {files}")
        
        for filename in files:
            match = re.match(pattern, filename)
            if match:
                condition = match.group(1)
                
                # 检查是否存在对应的注释文件
                anno_file = os.path.join(input_dir, f"{condition}_pasAnno.csv")
                anno_file_tsv = os.path.join(input_dir, f"{condition}_pasAnno.tsv")
                
                if os.path.exists(anno_file) or os.path.exists(anno_file_tsv):
                    conditions.append(condition)
                    logger.debug(f"找到条件: {condition}")
                else:
                    logger.warning(f"找到表达矩阵 {filename} 但没有对应的注释文件")
    except Exception as e:
        logger.error(f"获取处理条件时出错: {str(e)}")
        if logger.level <= logging.DEBUG:
            import traceback
            logger.debug(f"错误详细信息: {traceback.format_exc()}")
        raise
    
    return conditions

def main():
    """
    主函数
    """
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="Enhanced DECTDG-Analyzer: 增强版差异表达串联重复基因簇分析工具")
    
    # 输入文件
    parser.add_argument("--input_dir", required=True, help="输入文件目录")
    parser.add_argument("--gff_file", required=True, help="GFF文件路径")
    parser.add_argument("--gene_file", required=True, help="基因簇文件路径 (_processed_raw_genes_result)")
    parser.add_argument("--number_file", required=True, help="基因簇大小文件路径 (_processed_raw_numbers_result)")
    
    # 输出文件
    parser.add_argument("--output_dir", required=True, help="输出文件目录")
    
    # 参数设置
    parser.add_argument("--log2fc_threshold", type=float, default=1.0, help="Log2倍变阈值")
    parser.add_argument("--padj_threshold", type=float, default=0.05, help="校正P值阈值")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="使用的线程数")
    parser.add_argument("--verbose", action="store_true", help="显示详细日志")
    parser.add_argument("--debug", action="store_true", help="启用调试模式，记录更多信息")
    
    # 解析参数
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 设置日志文件
    log_file = os.path.join(args.output_dir, "enhanced_dectdg_analyzer.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    # 根据verbose参数设置日志级别
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.INFO)
    
    logger.info(f"Enhanced DECTDG-Analyzer v1.0 开始运行")
    logger.info(f"命令行参数: {sys.argv}")

    
    try:
        # 1. 加载基因长度信息
        gene_lengths = parse_gff(args.gff_file)
        
        # 2. 加载基因簇信息
        gene_clusters, cluster_numbers = load_gene_clusters(args.gene_file, args.number_file)
        
        # 3. 获取所有处理条件
        conditions = get_all_conditions(args.input_dir)
        if not conditions:
            logger.error("未找到任何处理条件")
            sys.exit(1)
        
        logger.info(f"找到 {len(conditions)} 个处理条件: {', '.join(conditions)}")
        
        # 4. 获取所有基因簇ID
        all_clusters = cluster_numbers['cluster_id'].unique().tolist()
        logger.info(f"共找到 {len(all_clusters)} 个基因簇")
        
        # 5. 多线程处理每个条件
        condition_results = []
        
        # 如果启用了调试模式，则不使用多线程
        if args.debug:
            logger.info("调试模式已启用，禁用多线程处理")
            for condition in conditions:
                try:
                    counts_file = os.path.join(args.input_dir, f"{condition}_counts_matrix.tsv")
                    anno_file = os.path.join(args.input_dir, f"{condition}_pasAnno.csv")
                    
                    # 先检查注释文件是否存在，不存在则尝试tsv格式
                    if not os.path.exists(anno_file):
                        alt_anno_file = os.path.join(args.input_dir, f"{condition}_pasAnno.tsv")
                        if os.path.exists(alt_anno_file):
                            anno_file = alt_anno_file
                    
                    # 加载数据
                    counts, anno = load_expression_data(counts_file, anno_file)
                    
                    # 计算TPM
                    tpm = calculate_tpm(counts, gene_lengths)
                    
                    # 获取处理组样本
                    treated_samples = anno[anno['condition'] == 'treated'].index.tolist()
                    
                    # 准备条件数据
                    condition_data = {
                        'counts': counts,
                        'anno': anno,
                        'gene_clusters': gene_clusters,
                        'cluster_numbers': cluster_numbers,
                        'log2fc_threshold': args.log2fc_threshold,
                        'padj_threshold': args.padj_threshold,
                        'tpm': tpm,
                        'condition_name': condition,
                        'threads_per_condition': 1,
                        'treated_samples': treated_samples
                    }
                    
                    # 处理条件
                    result = process_condition(condition_data)
                    condition_results.append(result)
                    logger.info(f"完成条件 {condition} 的处理")
                except Exception as e:
                    logger.error(f"处理条件 {condition} 时出错: {str(e)}")
                    if args.debug:
                        import traceback
                        logger.debug(f"错误详细信息: {traceback.format_exc()}")
        else:
            # 使用多线程处理
            threads_per_condition = max(1, min(2, args.threads // len(conditions)))
            logger.info(f"使用多线程处理，每个条件分配 {threads_per_condition} 个线程")
            
            with ProcessPoolExecutor(max_workers=args.threads) as executor:
                future_to_condition = {}
                
                for condition in conditions:
                    counts_file = os.path.join(args.input_dir, f"{condition}_counts_matrix.tsv")
                    anno_file = os.path.join(args.input_dir, f"{condition}_pasAnno.csv")
                    
                    # 先检查注释文件是否存在，不存在则尝试tsv格式
                    if not os.path.exists(anno_file):
                        alt_anno_file = os.path.join(args.input_dir, f"{condition}_pasAnno.tsv")
                        if os.path.exists(alt_anno_file):
                            anno_file = alt_anno_file
                    
                    # 加载数据
                    counts, anno = load_expression_data(counts_file, anno_file)
                    
                    # 计算TPM
                    tpm = calculate_tpm(counts, gene_lengths)
                    
                    # 获取处理组样本
                    treated_samples = anno[anno['condition'] == 'treated'].index.tolist()
                    
                    # 准备条件数据
                    condition_data = {
                        'counts': counts,
                        'anno': anno,
                        'gene_clusters': gene_clusters,
                        'cluster_numbers': cluster_numbers,
                        'log2fc_threshold': args.log2fc_threshold,
                        'padj_threshold': args.padj_threshold,
                        'tpm': tpm,
                        'condition_name': condition,
                        'threads_per_condition': threads_per_condition,
                        'treated_samples': treated_samples
                    }
                    
                    # 提交任务
                    future = executor.submit(process_condition, condition_data)
                    future_to_condition[future] = condition
                
                # 收集结果
                for future in as_completed(future_to_condition):
                    condition = future_to_condition[future]
                    try:
                        result = future.result()
                        condition_results.append(result)
                        logger.info(f"完成条件 {condition} 的处理")
                    except Exception as e:
                        logger.error(f"处理条件 {condition} 时出错: {str(e)}")
                        if args.debug:
                            import traceback
                            logger.debug(f"错误详细信息: {traceback.format_exc()}")
        
        # 检查是否有成功处理的条件
        if not condition_results:
            logger.error("没有成功处理任何条件，无法创建表达矩阵")
            sys.exit(1)
        
        # 6. 分析成员数量为2的DE-CTDGs中两个成员的表达差异
        member_expression_results = analyze_member_expression(condition_results)
        
        # 保存成员表达分析结果
        for condition_name, member_df in member_expression_results.items():
            member_file = os.path.join(args.output_dir, f"{condition_name}_member_expression.csv")
            member_df.to_csv(member_file, index=False)
            logger.info(f"{condition_name} 条件下成员数量为2的DE-CTDGs成员表达分析结果已保存到: {member_file}")
            
        # 7. 创建表达矩阵
        binary_matrix, tpm_matrix = create_expression_matrices(
            all_clusters, 
            condition_results, 
            cluster_numbers
        )
        
        # 8. 保存结果
        binary_output = os.path.join(args.output_dir, "binary_expression_matrix.csv")
        tpm_output = os.path.join(args.output_dir, "tpm_expression_matrix.csv")
        
        binary_matrix.to_csv(binary_output)
        tpm_matrix.to_csv(tpm_output)
        
        logger.info(f"01分布矩阵已保存到: {binary_output}")
        logger.info(f"TPM表达量矩阵已保存到: {tpm_output}")
        
        # 9. 为每个条件保存差异表达基因簇列表
        for result in condition_results:
            condition_name = result['condition_name']
            deg_clusters = result['deg_clusters']
            
            output_file = os.path.join(args.output_dir, f"{condition_name}_DEG_clusters.txt")
            with open(output_file, 'w') as f:
                for cluster in deg_clusters:
                    f.write(f"{cluster}\n")
            
            logger.info(f"{condition_name} 条件下的差异表达基因簇列表已保存到: {output_file}")
        
        # 10. 为每个条件保存详细的差异表达结果
        for result in condition_results:
            condition_name = result['condition_name']
            results_cluster = result['results_cluster']
            
            # 保存详细结果
            detailed_output = os.path.join(args.output_dir, f"{condition_name}_detailed_results.csv")
            results_cluster.to_csv(detailed_output)
            logger.info(f"{condition_name} 条件下的详细差异表达结果已保存到: {detailed_output}")
        
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"Enhanced DECTDG-Analyzer运行完成，总用时: {execution_time:.2f} 秒")
    
    except Exception as e:
        logger.error(f"程序执行过程中出错: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()