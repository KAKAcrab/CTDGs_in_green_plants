#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#########################################################################################################################
# 脚本名称: GO_enrichment_goatools_fixed.py
# 功能描述: 使用GOATOOLS对多物种不同逆境条件下的DE-CTDGs进行GO富集分析，与R版本topGO结果进行交叉验证
# 版本: v1.0
# 日期: 2025-05-04
#########################################################################################################################

import os
import sys
import time
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
from functools import partial
import logging
import threading

# 导入GOATOOLS相关库
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_ncbi_gene2go

# 设置日志
def setup_logger(log_file=None):
    """设置日志功能"""
    logger = logging.getLogger('GO_enrichment')
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    # 控制台输出
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # 文件输出
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

logger = setup_logger()

#########################################################################################################################
# 辅助函数 (Helper functions)
#########################################################################################################################

def read_gene2go(file_path):
    """读取gene2go映射文件, 转换为GOATOOLS需要的格式"""
    logger.info(f"读取gene2go映射文件: {file_path}")
    
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        # 读取原始gene2go文件
        gene2go_data = {}
        with open(file_path, 'r') as f:
            for line in f:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                gene_id = parts[0]
                go_ids = parts[1].split(',')
                
                # 确保GO ID格式正确
                go_ids = [go_id.strip() for go_id in go_ids if go_id.strip()]
                
                if gene_id not in gene2go_data:
                    gene2go_data[gene_id] = set()
                
                gene2go_data[gene_id].update(go_ids)
        
        # 注意：这里保持为集合类型而不是转换为列表
        # GOATOOLS期望gene_to_go中的值是集合而不是列表
        gene_to_go = gene2go_data  # 直接使用集合类型
        
        logger.info(f"成功加载 {len(gene2go_data)} 个基因的GO注释")
        return gene2go_data, gene_to_go
    
    except Exception as e:
        logger.error(f"错误: 读取gene2go映射文件时出错: {str(e)}")
        return {}, {}

def read_all_genes(file_path):
    """读取所有基因列表"""
    logger.info(f"读取所有基因列表: {file_path}")
    
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        all_genes = []
        with open(file_path, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    all_genes.append(gene)
        
        logger.info(f"成功读取 {len(all_genes)} 个基因")
        return all_genes
    
    except Exception as e:
        logger.error(f"错误: 读取所有基因列表时出错: {str(e)}")
        return []

def read_ctdg_genes(file_path):
    """读取基因簇成员数据"""
    logger.info(f"读取CTDG基因数据: {file_path}")
    
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        # 检查第一行是否为标题行
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            has_header = "gene" in first_line.lower()
        
        if has_header:
            logger.info("检测到标题行，将跳过第一行")
            genes_data = pd.read_csv(file_path, sep='\t', header=None, skiprows=1)
        else:
            genes_data = pd.read_csv(file_path, sep='\t', header=None)
        
        # 添加列名
        genes_data.columns = ["gene_id", "species_id", "chromosome_id", "cluster_id", "position_in_cluster"]
        
        logger.info(f"成功读取 {len(genes_data)} 条基因簇成员数据")
        return genes_data
    
    except Exception as e:
        logger.error(f"错误: 读取CTDG基因数据时出错: {str(e)}")
        return pd.DataFrame()

def read_expression_matrix(file_path):
    """读取表达矩阵"""
    logger.info(f"读取表达矩阵: {file_path}")
    
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        # 读取表达矩阵
        expr_matrix = pd.read_csv(file_path)
        
        # 检查数据格式
        if expr_matrix.shape[1] < 2:
            raise ValueError("表达矩阵格式不正确，至少需要2列")
        
        logger.info(f"成功读取表达矩阵，包含 {expr_matrix.shape[0]} 个基因簇和 {expr_matrix.shape[1]-1} 个条件")
        return expr_matrix
    
    except Exception as e:
        logger.error(f"错误: 读取表达矩阵时出错: {str(e)}")
        return pd.DataFrame()

def group_dectdgs_by_condition(expr_matrix, genes_data):
    """按条件分组DE-CTDGs"""
    logger.info("按条件分组DE-CTDGs")
    
    try:
        # 检查表达矩阵是否为空
        if expr_matrix.empty or expr_matrix.shape[1] <= 1:
            raise ValueError("表达矩阵为空或只有一列")
        
        # 获取条件列名（除了第一列以外的所有列）
        first_col_name = expr_matrix.columns[0]
        condition_cols = [col for col in expr_matrix.columns if col not in [first_col_name, "size"]]
        
        # 如果没有找到条件列，终止处理
        if not condition_cols:
            raise ValueError("未在表达矩阵中找到任何条件列")
        
        logger.info(f"找到 {len(condition_cols)} 个条件: {', '.join(condition_cols)}")
        
        # 以条件为键创建字典
        condition_genes = {}
        
        # 对每个条件，提取DE-CTDGs以及它们的成员基因
        for condition in condition_cols:
            # 提取当前条件下差异表达的基因簇
            de_clusters = expr_matrix.loc[expr_matrix[condition] == 1, first_col_name].tolist()
            
            if de_clusters:
                # 提取这些基因簇的所有成员基因
                mask = genes_data['cluster_id'].isin(de_clusters)
                de_genes = genes_data.loc[mask, 'gene_id'].unique().tolist()
                
                condition_genes[condition] = de_genes
                logger.info(f"条件 {condition} 下有 {len(de_clusters)} 个DE-CTDGs，包含 {len(de_genes)} 个唯一基因")
            else:
                condition_genes[condition] = []
                logger.info(f"条件 {condition} 下没有DE-CTDGs")
        
        return condition_genes
    
    except Exception as e:
        logger.error(f"错误: 按条件分组DE-CTDGs时出错: {str(e)}")
        return {}

def perform_go_enrichment(condition_name, condition_genes, all_genes, gene_to_go, go_dag, ontology="BP"):
    """对指定条件的基因进行GO富集分析"""
    logger.info(f"对条件 {condition_name} 的基因进行 {ontology} 本体富集分析")
    
    # 检查条件是否为空
    if not condition_genes:
        logger.warning(f"警告: 条件 {condition_name} 下没有基因，跳过富集分析")
        return pd.DataFrame()
    
    try:
        # 准备研究基因组和总体基因组 - 确保是集合类型
        study_genes = set(condition_genes)
        population = set(all_genes)
        
        # 确保研究基因组是总体基因组的子集
        if not study_genes.issubset(population):
            logger.warning("警告: 研究基因集不是总体基因集的子集，正在调整...")
            study_genes = study_genes.intersection(population)
            if not study_genes:
                logger.warning(f"警告: 调整后研究基因集为空，跳过富集分析")
                return pd.DataFrame()
        
        # 创建合适格式的关联字典 - 确保格式正确
        associations = {}
        for gene, go_ids in gene_to_go.items():
            if gene in population:  # 只包含总体基因集中的基因
                associations[gene] = go_ids  # 保持为集合类型
        
        # 创建GOEnrichmentStudy对象 - 使用位置参数和正确的数据格式
        goeaobj = GOEnrichmentStudy(
            population,           # 背景基因集
            associations,         # 基因到GO的关联字典
            go_dag,               # GO DAG对象
            propagate_counts=True,
            alpha=0.05,
            methods=['fdr_bh']
        )
        
        # 过滤只分析指定的本体类型
        namespaces = [ontology]
        
        # 执行富集分析
        goea_results = goeaobj.run_study(study_genes, namespaces=namespaces)
        
        # 筛选显著结果
        significant_results = [r for r in goea_results if r.p_fdr_bh < 0.05]
        
        # 如果没有显著结果，返回空DataFrame
        if not significant_results:
            logger.warning(f"警告: 条件 {condition_name} 的 {ontology} 富集分析没有找到显著的GO条目")
            return pd.DataFrame()
        
        # 转换为DataFrame格式
        results_df = pd.DataFrame([
            {
                'GO.ID': r.GO,
                'Term': r.name,
                'Annotated': r.pop_count,  # 总体中注释的基因数
                'Significant': r.study_count,  # 研究集中注释的基因数
                'Expected': r.pop_count * len(study_genes) / len(population),  # 期望值
                'p_uncorrected': r.p_uncorrected,  # 未校正的p值
                'p_fdr_bh': r.p_fdr_bh,  # FDR校正的p值
                'Group': condition_name,
                'ontology': ontology
            }
            for r in significant_results
        ])
        
        # 按p值排序
        results_df = results_df.sort_values('p_fdr_bh')
        
        return results_df
    
    except Exception as e:
        logger.error(f"错误: 在进行条件 {condition_name} 的 {ontology} 富集分析时发生错误: {str(e)}")
        return pd.DataFrame()

def perform_enrichment_for_all_conditions(condition_genes, all_genes, gene_to_go, go_dag, ontologies=None):
    """串行执行所有条件的GO富集分析"""
    if ontologies is None:
        ontologies = ["BP", "MF", "CC"]
    
    logger.info(f"开始为所有条件执行富集分析，选择的本体类型: {', '.join(ontologies)}")
    
    all_results = pd.DataFrame()
    
    # 对每种本体类型依次处理
    for ontology in ontologies:
        logger.info(f"开始处理本体类型: {ontology}")
        
        # 准备处理所需的参数
        condition_names = list(condition_genes.keys())
        
        # 创建结果列表
        results = []
        
        # 串行执行富集分析，避免多进程递归错误
        for condition_name in condition_names:
            condition_gene_list = condition_genes[condition_name]
            result = perform_go_enrichment(
                condition_name, 
                condition_gene_list, 
                all_genes, 
                gene_to_go, 
                go_dag, 
                ontology
            )
            results.append(result)
        
        # 过滤空结果
        non_empty_results = [df for df in results if not df.empty]
        
        if non_empty_results:
            # 合并所有结果
            combined_results = pd.concat(non_empty_results, ignore_index=True)
            all_results = pd.concat([all_results, combined_results], ignore_index=True)
        else:
            logger.warning(f"警告: 没有找到任何显著富集的 {ontology} GO条目")
    
    if not all_results.empty:
        return all_results
    else:
        logger.warning("警告: 所有本体类型均未发现显著富集的GO条目")
        # 返回一个空的数据框但包含所有需要的列
        return pd.DataFrame(columns=[
            'GO.ID', 'Term', 'Annotated', 'Significant', 'Expected',
            'p_uncorrected', 'p_fdr_bh', 'Group', 'ontology'
        ])

def create_go_dotplot(enrichment_data, output_dir, ontology=None, top_n=15):
    """创建GO富集点图 - 为每种本体类型创建单独的点图"""
    # 如果指定了本体类型，则筛选数据
    if ontology:
        enrichment_data = enrichment_data[enrichment_data['ontology'] == ontology]
        plot_title = f"{ontology} Ontology Enrichment"
        output_file = os.path.join(output_dir, f"DECTDGs_enrichment_of_{ontology}.pdf")
    else:
        plot_title = "GO Enrichment Analysis"
        output_file = os.path.join(output_dir, "DECTDGs_enrichment_all_ontologies.pdf")
    
    logger.info(f"创建GO富集点图: {output_file}")
    
    # 检查数据是否为空
    if enrichment_data.empty:
        logger.warning("警告: 没有数据，创建空白点图")
        # 创建一个空的点图
        plt.figure(figsize=(10, 10))
        plt.text(0.5, 0.5, "No significant GO terms found", 
                 horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        plt.savefig(output_file)
        plt.close()
        return None
    
    # 期望的条件顺序 - 自动检测所有条件并排序
    all_conditions = sorted(enrichment_data['Group'].unique())
    
    # 为每个条件选择前N个条目
    top_results = pd.DataFrame()
    
    for condition in all_conditions:
        # 提取当前条件的数据
        condition_data = enrichment_data[enrichment_data['Group'] == condition]
        
        if not condition_data.empty:
            # 按p值排序并选择前N个条目
            condition_top = condition_data.sort_values('p_fdr_bh').head(top_n)
            top_results = pd.concat([top_results, condition_top], ignore_index=True)
    
    # 检查是否还有数据
    if top_results.empty:
        logger.warning("警告: 筛选后没有数据，创建空白点图")
        plt.figure(figsize=(10, 10))
        plt.text(0.5, 0.5, "No significant GO terms found", 
                 horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        plt.savefig(output_file)
        plt.close()
        return None
    
    # 处理重名GO条目的问题
    # 创建点图所需的数据副本，以便修改显示格式而不影响原始数据
    plot_data = top_results.copy()
    
    # 初始化绘图专用的显示标签列
    plot_data['DisplayTerm'] = plot_data['Term']
    
    # 检查每个条件下是否有重复的Term
    for condition in plot_data['Group'].unique():
        condition_rows = plot_data[plot_data['Group'] == condition].index
        if len(condition_rows) > 0:
            # 获取该条件下的所有Term
            condition_terms = plot_data.loc[condition_rows, 'Term'].values
            
            # 检查是否有重复的Term
            if len(condition_terms) != len(set(condition_terms)):
                # 处理每个重复的Term
                duplicate_terms = set([term for term in condition_terms if list(condition_terms).count(term) > 1])
                
                for term in duplicate_terms:
                    # 找到有重复Term的所有行
                    term_rows = condition_rows[plot_data.loc[condition_rows, 'Term'] == term]
                    
                    # 给这些行的DisplayTerm添加GO.ID后缀
                    for row in term_rows:
                        # 使用完整GO ID作为标识
                        full_id = plot_data.loc[row, 'GO.ID']
                        # 更新显示标签
                        plot_data.loc[row, 'DisplayTerm'] = f"{plot_data.loc[row, 'Term']} [{full_id}]"
    
    # 将条件转换为因子，确保按期望顺序显示
    plot_data['Group'] = pd.Categorical(plot_data['Group'], categories=all_conditions, ordered=True)
    
    # 调整DisplayTerm顺序
    plot_data = plot_data.sort_values(['ontology', 'DisplayTerm'])
    unique_terms = plot_data['DisplayTerm'].unique()
    plot_data['DisplayTerm'] = pd.Categorical(plot_data['DisplayTerm'], categories=unique_terms, ordered=True)
    
    # 设置绘图参数
    plt.figure(figsize=(8.26, 9))
    
    # 创建点图
    g = sns.scatterplot(
        data=plot_data,
        x='Group',
        y='DisplayTerm',
        size='Significant',  # 基因数量决定点的大小
        hue='p_fdr_bh',      # p值决定点的颜色
        palette='coolwarm_r', # 颜色映射反转，低p值为红色
        sizes=(20, 200),      # 点的大小范围
        legend='brief'
    )
    
    # 设置图表标题和样式
    plt.title(plot_title, fontsize=14, fontweight='bold')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=0, ha='center', fontsize=8, fontweight='bold')
    plt.yticks(fontsize=8)
    
    # 调整图例
    plt.legend(title='Gene count', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 添加颜色条为p值
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax = plt.gca()
    axins = inset_axes(ax, width="5%", height="30%", loc='lower right', bbox_to_anchor=(1.05, 0., 1, 1), 
                      bbox_transform=ax.transAxes, borderpad=0)
    
    import matplotlib as mpl
    sm = plt.cm.ScalarMappable(cmap='coolwarm_r', norm=plt.Normalize(
        plot_data['p_fdr_bh'].min(), min(0.05, plot_data['p_fdr_bh'].max())
    ))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=axins)
    cbar.set_label('FDR adjusted p-value')
    
    # 保存图表
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    logger.info(f"点图已保存到: {output_file}")
    return output_file

def save_to_excel(enrichment_data, output_file):
    """保存结果到Excel文件"""
    logger.info(f"保存结果到Excel文件: {output_file}")
    
    # 检查数据是否为空
    if enrichment_data.empty:
        logger.warning("警告: 没有数据可以保存")
        # 创建一个空的工作簿
        empty_df = pd.DataFrame({'Message': ['No significant GO terms found']})
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            empty_df.to_excel(writer, sheet_name='Info', index=False)
        return
    
    # 创建ExcelWriter对象
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # 添加汇总工作表
        # 计算每个条件下各本体类型的GO条目数量
        summary_data = (
            enrichment_data
            .groupby(['Group', 'ontology'])
            .size()
            .reset_index(name='Term_Count')
            .pivot(index='Group', columns='ontology', values='Term_Count')
            .fillna(0)
            .astype(int)
            .reset_index()
        )
        
        # 写入汇总数据
        summary_data.to_excel(writer, sheet_name='Summary', index=False)
        
        # 为每种本体类型创建单独的工作表
        ontologies = enrichment_data['ontology'].unique()
        
        for ontology in ontologies:
            # 添加本体类型工作表
            sheet_name = f"{ontology}_Results"
            
            # 筛选当前本体类型的数据
            ontology_data = enrichment_data[enrichment_data['ontology'] == ontology].copy()
            
            # 创建一个包含所有条件TOP15结果的数据框
            conditions = sorted(ontology_data['Group'].unique())
            
            ontology_top_results = pd.DataFrame()
            
            for condition in conditions:
                # 提取条件数据并按p值排序
                condition_data = ontology_data[ontology_data['Group'] == condition].copy()
                condition_data = condition_data.sort_values('p_fdr_bh')
                
                # 选择前15个条目
                if len(condition_data) > 15:
                    condition_data = condition_data.head(15)
                
                # 添加到结果数据框
                ontology_top_results = pd.concat([ontology_top_results, condition_data], ignore_index=True)
            
            # 按条件排序
            ontology_top_results['Group'] = pd.Categorical(
                ontology_top_results['Group'], 
                categories=sorted(ontology_top_results['Group'].unique()),
                ordered=True
            )
            ontology_top_results = ontology_top_results.sort_values(['Group', 'p_fdr_bh'])
            
            # 写入数据
            ontology_top_results.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # 添加包含所有结果的工作表
        # 按条件和本体类型排序
        all_results = enrichment_data.copy()
        
        # 设置排序优先级
        all_results['Group'] = pd.Categorical(
            all_results['Group'],
            categories=sorted(all_results['Group'].unique()),
            ordered=True
        )
        all_results['ontology'] = pd.Categorical(
            all_results['ontology'],
            categories=['MF', 'CC', 'BP'],
            ordered=True
        )
        
        all_results = all_results.sort_values(['Group', 'ontology', 'p_fdr_bh'])
        
        # 写入所有数据
        all_results.to_excel(writer, sheet_name='All_Results', index=False)
    
    logger.info(f"结果已保存到: {output_file}")

def process_species(species_prefix, gene2go_file, all_genes_file, genes_file, expression_file, 
                    output_dir, go_obo_file, ontologies=None):
    """处理单个物种"""
    if ontologies is None:
        ontologies = ["BP", "MF", "CC"]
    
    logger.info(f"开始处理物种: {species_prefix}")
    
    # 创建物种输出目录
    species_output_dir = os.path.join(output_dir, f"{species_prefix}_GO_enrichment")
    os.makedirs(species_output_dir, exist_ok=True)
    
    # 设置输出文件路径
    output_xlsx = os.path.join(species_output_dir, f"{species_prefix}_DECTDGs_GO_enrichment.xlsx")
    
    # 读取GO本体
    logger.info(f"读取GO本体文件: {go_obo_file}")
    try:
        go_dag = GODag(go_obo_file)
        logger.info(f"成功加载GO本体，包含 {len(go_dag)} 个GO条目")
    except Exception as e:
        logger.error(f"错误: 读取GO本体时出错: {str(e)}")
        return False
    
    # 读取输入文件
    logger.info("读取输入文件")
    
    # 读取基因到GO映射文件
    gene2go_data, gene_to_go = read_gene2go(gene2go_file)
    if not gene2go_data:
        logger.error("错误: 无法读取gene2go映射文件，跳过物种处理")
        return False
    
    # 读取所有基因列表
    all_genes = read_all_genes(all_genes_file)
    if not all_genes:
        logger.error("错误: 无法读取所有基因列表，跳过物种处理")
        return False
    
    # 读取基因簇成员数据
    genes_data = read_ctdg_genes(genes_file)
    if genes_data.empty:
        logger.error("错误: 无法读取基因簇成员数据，跳过物种处理")
        return False
    
    # 读取01分布矩阵
    binary_matrix = read_expression_matrix(expression_file)
    if binary_matrix.empty:
        logger.error("错误: 无法读取01分布矩阵，跳过物种处理")
        return False
    
    # 按条件分组DE-CTDGs
    condition_genes = group_dectdgs_by_condition(binary_matrix, genes_data)
    if not condition_genes:
        logger.error("错误: 没有找到任何条件下的DE-CTDGs，跳过物种处理")
        return False
    
    # 执行GO富集分析，选择的本体类型
    enrichment_results = perform_enrichment_for_all_conditions(
        condition_genes, all_genes, gene_to_go, go_dag, ontologies
    )
    
    # 检查是否有结果
    if not enrichment_results.empty:
        logger.info("富集分析完成，开始创建可视化和保存结果")
        
        # 创建总体结果的点图
        create_go_dotplot(enrichment_results, species_output_dir, None, 15)
        
        # 为每种本体类型创建单独的点图
        for ontology in ontologies:
            create_go_dotplot(enrichment_results, species_output_dir, ontology, 15)
        
        # 保存结果到Excel文件
        save_to_excel(enrichment_results, output_xlsx)
        
        logger.info(f"处理完成，结果已保存到: {species_output_dir}")
        return True
    else:
        logger.warning("警告: 没有找到任何显著富集的GO条目")
        
        # 创建空的Excel文件
        empty_df = pd.DataFrame({'Message': ['No significant GO terms found']})
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            empty_df.to_excel(writer, sheet_name='Info', index=False)
        
        # 创建空的点图
        plt.figure(figsize=(10, 10))
        plt.text(0.5, 0.5, "No significant GO terms found", 
                 horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        
        output_file = os.path.join(species_output_dir, "DECTDGs_enrichment_all_ontologies.pdf")
        plt.savefig(output_file)
        plt.close()
        
        # 为每种本体类型创建空白点图
        for ontology in ontologies:
            plt.figure(figsize=(8.26, 9))
            plt.text(0.5, 0.5, f"No significant {ontology} GO terms found", 
                     horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
            output_file = os.path.join(species_output_dir, f"DECTDGs_enrichment_of_{ontology}.pdf")
            plt.savefig(output_file, dpi=300)
            plt.close()
        
        logger.info(f"处理完成，但没有找到显著结果，已创建空白输出文件: {species_output_dir}")
        return False

#########################################################################################################################
# 主函数 (Main function)
#########################################################################################################################

def main():
    """主函数"""
    logger.info("开始主函数执行")
    
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='GO富集分析工具，使用GOATOOLS')
    
    # 添加参数
    parser.add_argument('--species', type=str, default='Ppa', help='物种前缀，用于自动构建文件名')
    parser.add_argument('--gene2go', type=str, help='基因到GO映射文件路径')
    parser.add_argument('--all_genes', type=str, help='所有基因列表文件路径')
    parser.add_argument('--genes', type=str, help='基因簇成员数据文件路径')
    parser.add_argument('--expression', type=str, help='表达矩阵文件路径')
    parser.add_argument('--output', type=str, default='GO_enrichment_results', help='输出目录')
    parser.add_argument('--ontologies', type=str, default='BP,MF,CC', help='要分析的GO本体类型，用逗号分隔')
    parser.add_argument('--go_obo', type=str, default='go.obo', help='GO本体文件路径')
    parser.add_argument('--multi', action='store_true', help='是否处理多物种')
    parser.add_argument('--species_list', type=str, help='要处理的物种列表，用逗号分隔，仅在multi模式下有效')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置默认文件路径
    species_prefix = args.species
    gene2go_file = args.gene2go or f"{species_prefix}_gene2go"
    all_genes_file = args.all_genes or f"{species_prefix}_all_genes_list"
    genes_file = args.genes or f"{species_prefix}_processed_raw_genes_result"
    expression_file = args.expression or f"{species_prefix}_DE_CTDGs_matrix.csv"
    output_dir = args.output
    go_obo_file = args.go_obo
    
    # 解析本体类型
    ontologies = [ont.strip().upper() for ont in args.ontologies.split(',') if ont.strip()]
    # 验证本体类型
    valid_ontologies = ['BP', 'MF', 'CC']
    ontologies = [ont for ont in ontologies if ont in valid_ontologies]
    if not ontologies:
        logger.warning("警告: 未指定有效的本体类型，使用默认值 BP,MF,CC")
        ontologies = ['BP', 'MF', 'CC']
    
    # 创建总输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 检查是否处理多个物种
    if args.multi:
        # 处理多物种模式
        logger.info("多物种模式已激活")
        
        # 解析物种列表
        species_list = []
        if args.species_list:
            species_list = [sp.strip() for sp in args.species_list.split(',') if sp.strip()]
        else:
            # 尝试获取当前目录中所有符合条件的物种
            import glob
            matrix_files = glob.glob("*_DE_CTDGs_matrix.csv")
            species_list = [os.path.basename(f).replace('_DE_CTDGs_matrix.csv', '') for f in matrix_files]
        
        if not species_list:
            logger.error("错误: 未找到任何物种数据，终止处理")
            return
        
        logger.info(f"找到 {len(species_list)} 个物种: {', '.join(species_list)}")
        
        # 逐个处理每个物种
        success_count = 0
        for species in species_list:
            logger.info(f"开始处理物种: {species}")
            
            # 构建物种文件名
            sp_gene2go_file = f"{species}_gene2go"
            sp_all_genes_file = f"{species}_all_genes_list"
            sp_genes_file = f"{species}_processed_raw_genes_result"
            sp_expression_file = f"{species}_DE_CTDGs_matrix.csv"
            
            # 处理物种
            result = process_species(
                species, sp_gene2go_file, sp_all_genes_file, sp_genes_file,
                sp_expression_file, output_dir, go_obo_file, ontologies
            )
            
            if result:
                success_count += 1
        
        logger.info(f"多物种处理完成，成功处理 {success_count} 个物种，共 {len(species_list)} 个物种")
    else:
        # 单物种模式
        logger.info(f"分析物种: {species_prefix}，使用本体类型: {', '.join(ontologies)}")
        
        # 处理单个物种
        result = process_species(
            species_prefix, gene2go_file, all_genes_file, genes_file,
            expression_file, output_dir, go_obo_file, ontologies
        )
        
        if result:
            logger.info(f"物种 {species_prefix} 处理成功")
        else:
            logger.info(f"物种 {species_prefix} 处理完成，但可能没有发现显著富集")
    
    logger.info("脚本执行完成")

if __name__ == "__main__":
    main()