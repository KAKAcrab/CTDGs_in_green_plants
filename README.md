# CTDGs_in_green_plants

Comprehensive identification and analysis of clusters of tandemly duplicated genes reveals their contributions to adaptive evolution of green plants


**· Extract ctdgFinder results files: combines overlapping cluster results for all species listed in the files under CTDG_out directory and merges clusters with shared bases:**

```bash
python extract_result.py list ## get initial extracted results

## merge the "genes_result" of Aty, if not been merged completely, try 2-3 times until no "should have been merged" info showing
## round 1
python parse_result_of_CTDG_finder.py Aty
python validate_merge_result.py Aty
##round 2
mv Aty_processed_genes_results Aty_genes_result
python parse_result_of_CTDG_finder.py Aty
python validate_merge_result.py Aty
## round x
## ...

## we can also check out the partial overlapped information of the final merged genes_result
python overlap_analysis.py Aty
## finally we can generate the numbers_result file from the merged genes_file
python process_numbers_result.py Aty_processed_raw_genes_result Aty_processed_raw_numbers_result

```

**· Clustering rates analysis of different genes across 220 species**

Calculate the distribution matrix of clustering rates of genes of 7 traits in 220 species and obtain the protein sequences of all query
```bash
python analyze_plant_genes.py --result_dir ./ --clustered_dir ./ --species_order 220_spe_list_low2high_latin_abbr --gene_info 220_genes_parsed.csv --protein_file 220.pep --output_csv out_put/gene_distribution_stats.csv --output_seq_dir out_put
```

Select genes with the largest cluster change rate

```bash
Rscript ../custom_script/cluster_ratio_analyse.R -i ./gene_distribution_stats.csv -o ./cluster_ratio_analysis
```

Calculate the distribution matrix of clustering rates of different genes in different species and the bubble diagram of the first three traits

```bash
Rscript ../../custom_script/Bubble_plot_of_cluster_ratio.R -i gene_distribution_stats.csv -g group2spe.csv -n gene_family_name2_gene_ID.csv -o result
```





**· For each stress condition, caculate the DE-CTDGs according to the Mean Expression Difference and the Percentage of DEGs in Members, finally get the DE-CTDGs matrix for further analysis:**

```bash
## trimming
java -jar PATH/TO/trimmomatic-0.39.jar PE -threads 10 *.sra.fastq.gz *.sra.fastq.gz -baseout 1.clean_reads/* ILLUMINACLIP:PATH/TO/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

## quality control
fastqc -o . *.sra.fastq.gz

## align
hisat2 -t  -p 5 -x Aty -1 *.sra_1.fastq.gz -2 *.sra_2.fastq.gz -S *.sam
samtools view -S *.sam -b > *.bam
samtools sort -n *.bam -o *_nsorted.bam

## get reads count
htseq-count -n 10 -t gene -i ID -r name -f bam *_nsorted.bam Aty.gff3 > *.count

## calculate each conditon one by one and get the DE-CTDGs lists
python identify_DECTDG.py --gene_file Aty_processed_raw_genes_result --number_file Aty_processed_raw_numbers_result --counts_file cold_counts_matrix.tsv --anno_file cold_pasAnno.csv --output_file cold_DEG_clusters.txt

## merge the DE-CTDGs lists under all conditions into the final matrix
python create_expression_matrix.py all_CTDGs_list_of_Aty Aty_processed_raw_numbers_result expression_matrix.csv cold_DEG_clusters.txt drought_DEG_clusters.txt heat_DEG_clusters.txt light_DEG_clusters.txt salt_DEG_clusters.txt
```
