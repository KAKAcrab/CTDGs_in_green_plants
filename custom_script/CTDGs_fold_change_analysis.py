#!/usr/bin/env python3

"""
CTDG-FoldChange Analysis Script
Author: yyh
Date: 2025-05-20

计算不同大小real CTDGs和random CTDGs在不同条件下平均TPM的fold change

功能:
1. 自动读取输入目录下的所有相关文件
2. 自动检测并分析所有条件
3. 按成员数量(2, 3-5, 6-10, >10)将CTDGs分组
4. 计算每组真实CTDGs在处理前后的fold change (使用伪计数法，TPM+1)
5. 随机生成对应数量的随机CTDGs并计算fold change
6. 输出密度图和箱形图数据供下游使用R脚本可视化

使用方法:
  python ctdg_fold_change_analysis.py \
    --input_dir <input_directory> \
    --output_dir <output_directory> \
    [--random_count <count>]
"""

import os
import sys
import argparse
import random
import json
import csv
import glob
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

# 设置随机种子以保证结果可重复
random.seed(42)
np.random.seed(42)

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='分析CTDGs的fold change')
    
    parser.add_argument('--input_dir', required=True, help='包含所有输入文件的目录')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    parser.add_argument('--random_count', type=int, default=1000, 
                       help='每组随机CTDGs的数量 (默认: 1000)')
    
    return parser.parse_args()

def find_input_files(input_dir):
    """
    在输入目录中查找所有必要的输入文件
    
    参数:
        input_dir (str): 输入目录路径
        
    返回:
        tuple: (gff_file, genes_result, numbers_result, conditions_data)
            conditions_data是一个字典，键为条件名，值为(counts_matrix, sample_anno)
    """
    print(f"在目录 {input_dir} 中查找输入文件...")
    
    # 查找GFF文件
    gff_files = glob.glob(os.path.join(input_dir, "*.gff")) + glob.glob(os.path.join(input_dir, "*.gff3"))
    if not gff_files:
        raise FileNotFoundError(f"在 {input_dir} 中未找到GFF文件")
    gff_file = gff_files[0]
    print(f"找到GFF文件: {os.path.basename(gff_file)}")
    
    # 查找基因簇结果文件
    genes_result_files = glob.glob(os.path.join(input_dir, "*processed_raw_genes_result*"))
    if not genes_result_files:
        raise FileNotFoundError(f"在 {input_dir} 中未找到基因簇结果文件")
    genes_result = genes_result_files[0]
    print(f"找到基因簇结果文件: {os.path.basename(genes_result)}")
    
    # 查找基因簇数量结果文件
    numbers_result_files = glob.glob(os.path.join(input_dir, "*processed_raw_numbers_result*"))
    if not numbers_result_files:
        raise FileNotFoundError(f"在 {input_dir} 中未找到基因簇数量结果文件")
    numbers_result = numbers_result_files[0]
    print(f"找到基因簇数量结果文件: {os.path.basename(numbers_result)}")
    
    # 查找条件相关文件
    counts_matrices = glob.glob(os.path.join(input_dir, "*_counts_matrix.tsv"))
    sample_annos = glob.glob(os.path.join(input_dir, "*_pasAnno.csv"))
    
    # 提取条件名
    condition_pattern = re.compile(r'(.+)_counts_matrix\.tsv')
    conditions = set()
    for counts_file in counts_matrices:
        match = condition_pattern.search(os.path.basename(counts_file))
        if match:
            conditions.add(match.group(1))
    
    # 检查每个条件是否有对应的注释文件
    conditions_data = {}
    for condition in conditions:
        counts_file = os.path.join(input_dir, f"{condition}_counts_matrix.tsv")
        anno_file = os.path.join(input_dir, f"{condition}_pasAnno.csv")
        
        if os.path.exists(counts_file) and os.path.exists(anno_file):
            conditions_data[condition] = (counts_file, anno_file)
            print(f"找到条件 {condition} 的文件")
        else:
            print(f"警告: 条件 {condition} 缺少必要的文件")
    
    if not conditions_data:
        raise FileNotFoundError(f"在 {input_dir} 中未找到有效的条件数据")
    
    return gff_file, genes_result, numbers_result, conditions_data

def read_gff(gff_file):
    """
    读取GFF文件，提取基因信息
    
    参数:
        gff_file (str): GFF文件路径
        
    返回:
        dict: 基因ID -> (染色体ID, 起始位置, 结束位置)
    """
    print(f"读取GFF文件: {gff_file}")
    genes = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            if fields[2] == "gene":
                # 提取基因ID
                attr = fields[8]
                gene_id = None
                
                # 查找ID字段
                for item in attr.split(';'):
                    if item.startswith('ID='):
                        gene_id = item[3:].split(',')[0]
                        break
                
                if gene_id:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    genes[gene_id] = (chrom, start, end)
    
    print(f"成功读取 {len(genes)} 个基因")
    return genes

def read_genes_result(genes_result_file):
    """
    读取成簇基因信息文件
    
    参数:
        genes_result_file (str): 成簇基因信息文件路径
        
    返回:
        dict: 基因簇ID -> 基因ID列表
    """
    print(f"读取成簇基因信息: {genes_result_file}")
    ctdg_genes = defaultdict(list)
    
    with open(genes_result_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 5:
                continue
                
            gene_id = row[0]
            ctdg_id = row[3]
            
            # 将基因添加到对应的基因簇
            ctdg_genes[ctdg_id].append(gene_id)
    
    # 去除重复的基因ID
    for ctdg_id in ctdg_genes:
        ctdg_genes[ctdg_id] = list(dict.fromkeys(ctdg_genes[ctdg_id]))
    
    print(f"成功读取 {len(ctdg_genes)} 个基因簇")
    return ctdg_genes

def read_numbers_result(numbers_result_file):
    """
    读取基因簇大小信息文件
    
    参数:
        numbers_result_file (str): 基因簇大小信息文件路径
        
    返回:
        dict: 基因簇ID -> 成员数量
    """
    print(f"读取基因簇大小信息: {numbers_result_file}")
    ctdg_sizes = {}
    
    with open(numbers_result_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 4:
                continue
                
            ctdg_id = row[2]
            
            # 确保成员数量是整数
            try:
                size = int(row[3])
                ctdg_sizes[ctdg_id] = size
            except ValueError:
                # 跳过无法转换为整数的值
                continue
    
    print(f"成功读取 {len(ctdg_sizes)} 个基因簇大小信息")
    return ctdg_sizes

def read_expression_data(counts_matrix_file, sample_anno_file):
    """
    读取表达量矩阵和样本注释文件
    
    参数:
        counts_matrix_file (str): 表达量矩阵文件路径
        sample_anno_file (str): 样本注释文件路径
        
    返回:
        tuple: (表达矩阵DataFrame, 样本注释DataFrame)
    """
    print(f"读取表达量矩阵: {os.path.basename(counts_matrix_file)}")
    # 读取表达量矩阵
    expr_df = pd.read_csv(counts_matrix_file, sep='\t', index_col=0)
    
    print(f"读取样本注释: {os.path.basename(sample_anno_file)}")
    # 读取样本注释
    anno_df = pd.read_csv(sample_anno_file)
    
    # 确保样本注释中的样本名与表达矩阵匹配
    common_samples = set(anno_df['sample']) & set(expr_df.index)
    if len(common_samples) == 0:
        raise ValueError(f"样本注释和表达矩阵中没有共同的样本: {sample_anno_file}")
    
    # 过滤表达矩阵中不在注释中的样本
    expr_df = expr_df.loc[expr_df.index.isin(common_samples)]
    
    print(f"成功读取 {expr_df.shape[1]} 个基因的表达量，包含 {expr_df.shape[0]} 个样本")
    return expr_df, anno_df

def group_ctdgs_by_size(ctdg_genes, ctdg_sizes):
    """
    按大小分组CTDGs
    
    参数:
        ctdg_genes (dict): 基因簇ID -> 基因ID列表
        ctdg_sizes (dict): 基因簇ID -> 成员数量
        
    返回:
        dict: 大小组 -> 基因簇ID列表
    """
    print("按大小分组CTDGs")
    size_groups = {
        "member=2": [],
        "3<member<5": [],
        "6<member<10": [],
        "member>10": []
    }
    
    for ctdg_id, genes in ctdg_genes.items():
        # 使用实际的基因列表长度而不是ctdg_sizes中的值
        size = len(genes)
        
        # 分配到对应的大小组
        if size == 2:
            size_groups["member=2"].append(ctdg_id)
        elif 3 <= size <= 5:
            size_groups["3<member<5"].append(ctdg_id)
        elif 6 <= size <= 10:
            size_groups["6<member<10"].append(ctdg_id)
        elif size > 10:
            size_groups["member>10"].append(ctdg_id)
    
    # 打印每组的数量
    for group, ctdgs in size_groups.items():
        print(f"  {group}: {len(ctdgs)} 个基因簇")
    
    return size_groups

def calculate_average_tpm(gene_ids, expr_df, condition_samples):
    """
    计算基因簇的平均表达量
    
    参数:
        gene_ids (list): 基因ID列表
        expr_df (DataFrame): 表达量矩阵
        condition_samples (list): 特定条件下的样本列表
        
    返回:
        float: 基因簇的平均TPM值（如果没有表达数据，返回0）
    """
    # 过滤存在于表达矩阵中的基因
    valid_genes = [g for g in gene_ids if g in expr_df.columns]
    
    if not valid_genes or not condition_samples:
        return 0.0
    
    # 确保所有样本都在表达矩阵中
    valid_samples = [s for s in condition_samples if s in expr_df.index]
    
    if not valid_samples:
        return 0.0
    
    # 提取表达量数据
    gene_expr = expr_df.loc[valid_samples, valid_genes]
    
    # 计算每个样本的平均表达量
    avg_expr = gene_expr.mean(axis=1)
    
    # 计算所有样本的平均值，如果有NaN值，则返回0
    mean_expr = avg_expr.mean()
    return 0.0 if pd.isna(mean_expr) else mean_expr

def generate_random_ctdgs(genes_info, size_range, count=1000):
    """
    生成随机CTDGs，每个CTDG的大小在指定范围内随机选择
    
    参数:
        genes_info (dict): 基因信息字典 (基因ID -> (染色体, 起始位置, 结束位置))
        size_range (tuple or list or int): 基因簇大小范围 
                                          (如(3,5)表示3到5之间，2表示固定为2)
        count (int): 生成的随机CTDGs数量
        
    返回:
        list: 随机CTDGs列表，每个元素是一个基因ID列表
    """
    random_ctdgs = []
    
    # 按染色体分组基因
    chrom_genes = defaultdict(list)
    for gene_id, (chrom, start, end) in genes_info.items():
        chrom_genes[chrom].append((gene_id, start, end))
    
    # 对每个染色体上的基因进行排序（按位置）
    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda x: x[1])
    
    # 生成随机CTDGs
    attempts = 0
    max_attempts = count * 10  # 最大尝试次数
    
    while len(random_ctdgs) < count and attempts < max_attempts:
        attempts += 1
        
        # 为当前CTDG随机选择一个大小
        if isinstance(size_range, (tuple, list)) and len(size_range) == 2:
            # 范围内随机选择
            ctdg_size = random.randint(size_range[0], size_range[1])
        elif isinstance(size_range, int):
            # 固定大小
            ctdg_size = size_range
        else:
            # 默认大小（如果参数格式不正确）
            ctdg_size = 2
        
        # 随机选择一个染色体
        chroms = list(chrom_genes.keys())
        if not chroms:
            continue
            
        chrom = random.choice(chroms)
        genes = chrom_genes[chrom]
        
        # 确保染色体上有足够的基因
        if len(genes) < ctdg_size:
            continue
        
        # 随机选择起始位置
        start_idx = random.randint(0, len(genes) - ctdg_size)
        
        # 选择连续的基因
        random_ctdg = [genes[start_idx + i][0] for i in range(ctdg_size)]
        
        # 添加到结果列表
        random_ctdgs.append(random_ctdg)
    
    return random_ctdgs[:count]  # 确保返回正确数量的随机CTDGs

def process_condition(condition, counts_matrix_file, sample_anno_file, ctdg_genes, size_groups, genes_info, output_dir, random_count):
    """
    处理特定条件下的CTDGs fold change
    
    参数:
        condition (str): 条件名称
        counts_matrix_file (str): 表达量矩阵文件路径
        sample_anno_file (str): 样本注释文件路径
        ctdg_genes (dict): 基因簇ID -> 基因ID列表
        size_groups (dict): 大小组 -> 基因簇ID列表
        genes_info (dict): 基因信息
        output_dir (str): 输出目录
        random_count (int): 随机CTDGs的数量
    """
    print(f"\n处理条件: {condition}")
    
    # 创建条件输出目录
    condition_dir = os.path.join(output_dir, "visualization_data", condition)
    os.makedirs(condition_dir, exist_ok=True)
    
    # 读取表达量数据和样本注释
    try:
        expr_df, anno_df = read_expression_data(counts_matrix_file, sample_anno_file)
    except Exception as e:
        print(f"处理条件 {condition} 时出错: {str(e)}")
        return
    
    # 获取处理前后的样本
    treated_samples = anno_df[anno_df['condition'] == 'treated']['sample'].tolist()
    untreated_samples = anno_df[anno_df['condition'] == 'untreated']['sample'].tolist()
    
    print(f"  处理前样本数: {len(untreated_samples)}")
    print(f"  处理后样本数: {len(treated_samples)}")
    
    if not treated_samples or not untreated_samples:
        print(f"  警告: 条件 {condition} 缺少处理前或处理后的样本，跳过")
        return
    
    # 存储密度图和箱形图数据
    density_data = []
    boxplot_data = []
    
    # 按大小组处理
    for size_group, ctdg_ids in size_groups.items():
        if not ctdg_ids:
            print(f"  跳过空组 {size_group}")
            continue
            
        print(f"  处理大小组: {size_group} (包含 {len(ctdg_ids)} 个CTDGs)")
        
        # 计算真实CTDGs的fold change
        real_fold_changes = []
        
        for ctdg_id in ctdg_ids:
            genes = ctdg_genes[ctdg_id]
            
            # 计算处理前后的平均表达量
            treated_expr = calculate_average_tpm(genes, expr_df, treated_samples)
            untreated_expr = calculate_average_tpm(genes, expr_df, untreated_samples)
            
            # 应用伪计数（对处理前后的表达量都加1）
            treated_expr_pseudo = treated_expr + 1.0
            untreated_expr_pseudo = untreated_expr + 1.0
            
            # 计算fold change (log2)
            fold_change = np.log2(treated_expr_pseudo / untreated_expr_pseudo)
            real_fold_changes.append(fold_change)
        
        if not real_fold_changes:
            print(f"  警告: {size_group} 没有有效的fold change值，跳过")
            continue
            
        print(f"    计算得到 {len(real_fold_changes)} 个真实CTDGs的fold change")
        
        # 添加真实CTDGs的密度图数据
        for fc in real_fold_changes:
            density_data.append({
                "size_group": size_group,
                "group": "CTDG",
                "log2_fold_change": fc
            })
            
            boxplot_data.append({
                "size_group": size_group,
                "group": "CTDG",
                "log2_fold_change": fc
            })
        
        # 设置随机CTDGs的大小范围
        if size_group == "member=2":
            size_range = 2  # 固定为2
        elif size_group == "3<member<5":
            size_range = (3, 5)  # 3到5之间随机
        elif size_group == "6<member<10":
            size_range = (6, 10)  # 6到10之间随机
        else:  # member>10
            size_range = (11, 20)  # 11到20之间随机（限制最大随机大小为20）
        
        # 生成随机CTDGs
        random_ctdgs = generate_random_ctdgs(genes_info, size_range, count=random_count)
        print(f"    生成 {len(random_ctdgs)} 个随机CTDGs (大小范围: {size_range})")
        
        # 计算随机CTDGs的fold change
        random_fold_changes = []
        
        for genes in random_ctdgs:
            # 计算处理前后的平均表达量
            treated_expr = calculate_average_tpm(genes, expr_df, treated_samples)
            untreated_expr = calculate_average_tpm(genes, expr_df, untreated_samples)
            
            # 应用伪计数（对处理前后的表达量都加1）
            treated_expr_pseudo = treated_expr + 1.0
            untreated_expr_pseudo = untreated_expr + 1.0
            
            # 计算fold change (log2)
            fold_change = np.log2(treated_expr_pseudo / untreated_expr_pseudo)
            random_fold_changes.append(fold_change)
        
        print(f"    计算得到 {len(random_fold_changes)} 个随机CTDGs的fold change")
        
        # 添加随机CTDGs的密度图数据
        for fc in random_fold_changes:
            density_data.append({
                "size_group": size_group,
                "group": "Random",
                "log2_fold_change": fc
            })
            
            boxplot_data.append({
                "size_group": size_group,
                "group": "Random",
                "log2_fold_change": fc
            })
    
    # 输出密度图数据
    if density_data:
        density_file = os.path.join(condition_dir, "density_data.csv")
        pd.DataFrame(density_data).to_csv(density_file, index=False)
        print(f"  密度图数据已保存至: {density_file}")
    
    # 输出箱形图数据
    if boxplot_data:
        boxplot_file = os.path.join(condition_dir, "boxplot_data.csv")
        pd.DataFrame(boxplot_data).to_csv(boxplot_file, index=False)
        print(f"  箱形图数据已保存至: {boxplot_file}")

def write_condition_list(conditions, output_dir):
    """
    将检测到的条件列表写入文件
    
    参数:
        conditions (list): 条件列表
        output_dir (str): 输出目录
    """
    conditions_file = os.path.join(output_dir, "conditions.txt")
    with open(conditions_file, 'w') as f:
        for condition in conditions:
            f.write(f"{condition}\n")
    print(f"条件列表已保存至: {conditions_file}")

def main():
    """主函数"""
    args = parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    vis_dir = os.path.join(args.output_dir, "visualization_data")
    os.makedirs(vis_dir, exist_ok=True)
    
    # 查找输入文件
    try:
        gff_file, genes_result, numbers_result, conditions_data = find_input_files(args.input_dir)
    except Exception as e:
        print(f"查找输入文件时出错: {str(e)}")
        sys.exit(1)
    
    # 读取基本数据
    genes_info = read_gff(gff_file)
    ctdg_genes = read_genes_result(genes_result)
    ctdg_sizes = read_numbers_result(numbers_result)
    
    # 按大小分组CTDGs
    size_groups = group_ctdgs_by_size(ctdg_genes, ctdg_sizes)
    
    # 将条件列表写入文件
    write_condition_list(list(conditions_data.keys()), args.output_dir)
    
    # 处理每个条件
    for condition, (counts_matrix_file, sample_anno_file) in conditions_data.items():
        process_condition(condition, counts_matrix_file, sample_anno_file, 
                        ctdg_genes, size_groups, genes_info, args.output_dir, args.random_count)
    
    print("\n分析完成!")

if __name__ == '__main__':
    main()