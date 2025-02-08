import pandas as pd
import argparse
import os

def load_deg_clusters(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f)

def load_cluster_sizes(file_path):
    sizes = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                cluster_id = parts[2]
                size = int(parts[3])
                if size ==2:
                    sizes[cluster_id] = "member=2";
                elif 2 < size <= 5:
                    sizes[cluster_id] = "3<member<5"
                elif 5 < size <= 10:
                    sizes[cluster_id] = "6<member<10"
                else :
                    sizes[cluster_id] = "member>10"

    return sizes

def create_expression_matrix(all_clusters_file, treatment_files, cluster_sizes_file, output_file):
    # 读取所有基因簇ID
    with open(all_clusters_file, 'r') as f:
        all_clusters = [line.strip() for line in f]

    # 读取基因簇大小信息
    cluster_sizes = load_cluster_sizes(cluster_sizes_file)

    # 创建一个DataFrame，初始值都为0
    columns = [os.path.basename(f).split('_')[0] for f in treatment_files]
    matrix = pd.DataFrame(0, index=all_clusters, columns=columns)

    # 添加size列
    matrix['size'] = matrix.index.map(lambda x: cluster_sizes.get(x, 0))

    # 读取每个处理条件下的差异表达基因簇，并在矩阵中标记为1
    for file in treatment_files:
        treatment = os.path.basename(file).split('_')[0]
        deg_clusters = load_deg_clusters(file)
        matrix.loc[matrix.index.isin(deg_clusters), treatment] = 1

    # 将cluster列设为索引
    matrix.index.name = 'cluster'

    # 重新排列列，使得'size'列在最前面
    cols = matrix.columns.tolist()
    cols = ['size'] + [col for col in cols if col != 'size']
    matrix = matrix[cols]

    # 保存结果
    matrix.to_csv(output_file)
    print(f"Expression matrix has been saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Create expression matrix for differentially expressed gene clusters")
    parser.add_argument("all_clusters_file", help="File containing all cluster IDs")
    parser.add_argument("cluster_sizes_file", help="File containing cluster sizes (Aty_processed_raw_numbers_result)")
    parser.add_argument("output_file", help="Output file path for the expression matrix")
    parser.add_argument("treatment_files", nargs='+', help="Files containing differentially expressed clusters for each treatment")

    args = parser.parse_args()

    create_expression_matrix(args.all_clusters_file, args.treatment_files, args.cluster_sizes_file, args.output_file)

if __name__ == "__main__":
    main()
