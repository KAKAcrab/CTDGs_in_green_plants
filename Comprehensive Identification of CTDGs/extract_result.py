#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
from typing import List
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

class ResultProcessor:
    def __init__(self, spe_list_file: str, num_threads: int = 4):
        """初始化处理器"""
        self.work_path = Path.cwd()
        self.species = self._read_species_list(spe_list_file)
        self.num_threads = num_threads
        self.write_lock = Lock()

    def _read_species_list(self, spe_list_file: str) -> List[str]:
        """读取物种列表"""
        try:
            with open(spe_list_file, 'r') as f:
                species = [line.strip() for line in f if line.strip()]
            print(f"已读取物种列表:\n{', '.join(species)}")
            return species
        except FileNotFoundError:
            print(f"错误：找不到物种列表文件 {spe_list_file}")
            exit(1)

    def _process_file(self, file_path: Path, family_name: str, file_type: str) -> List[str]:
        """处理单个文件"""
        results = []
        try:
            with open(file_path, "r") as infile:
                next(infile)  # 跳过标题行
                for line in infile:
                    fields = line.strip().split(",")
                    if file_type == "genes" and len(fields) >= 5:
                        acc, spe, chrom, cluster, order = fields[0], fields[1], fields[2], fields[3], fields[-1]
                        results.append(f"{acc}\t{spe}\t{chrom}\t{cluster}\t{order}")
                    elif file_type == "numbers" and len(fields) >= 8:
                        spe, chrom, cluster, gene_dup, start, end, length, p_95 = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]
                        results.append(f"{spe}\t{family_name}\t{chrom}\t{cluster}\t{gene_dup}\t{start}\t{end}\t{length}\t{p_95}")
        except Exception as e:
            print(f"处理文件 {file_path} 时出错：{str(e)}")
        return results

    def _process_species_directory(self, species: str, file_type: str) -> List[str]:
        """处理单个物种目录"""
        results = []
        species_path = self.work_path / species
        if not species_path.is_dir():
            print(f"警告：找不到物种目录 {species}")
            return results

        for family in species_path.iterdir():
            if not family.is_dir():
                continue

            report_path = family / "report"
            file_suffix = f"{'genes' if file_type == 'genes' else 'numbers'}_clean.csv"
            file_name = f"{family.name}_{file_suffix}"
            file_path = report_path / file_name

            if not file_path.exists():
                continue

            results.extend(self._process_file(file_path, family.name, file_type))

        return results

    def process_results(self, file_type: str, output_file: str):
        """并行处理结果文件"""
        print(f"正在处理{file_type}结果...")
        
        with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
            future_to_species = {
                executor.submit(self._process_species_directory, species, file_type): species
                for species in self.species
            }

            # 使用write_lock确保写入操作的线程安全
            with open(output_file, "w") as out:
                for future in as_completed(future_to_species):
                    species = future_to_species[future]
                    try:
                        results = future.result()
                        with self.write_lock:
                            for result in results:
                                out.write(f"{result}\n")
                    except Exception as e:
                        print(f"处理物种 {species} 时出错：{str(e)}")

        print(f"{file_type}结果已保存至：{output_file}")

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='处理CTDG输出结果')
    parser.add_argument('spe_list', help='物种列表文件路径')
    parser.add_argument('--genes_out', default='genes_result.txt', help='genes结果输出文件路径')
    parser.add_argument('--numbers_out', default='numbers_result.txt', help='numbers结果输出文件路径')
    parser.add_argument('--threads', type=int, default=4, help='并行处理的线程数')
    return parser.parse_args()

def main():
    """主函数"""
    args = parse_args()
    
    try:
        processor = ResultProcessor(args.spe_list, args.threads)
        
        # 处理genes结果
        processor.process_results("genes", args.genes_out)
        
        # 处理numbers结果
        processor.process_results("numbers", args.numbers_out)
        
        print("处理完成！")
    except Exception as e:
        print(f"程序执行出错：{str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
