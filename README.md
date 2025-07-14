# CTDGs_in_green_plants

Comprehensive identification and analysis of clusters of tandemly duplicated genes (CTDGs) reveals their contributions to adaptive evolution of green plants. Input files for ctdgFinder of 220 plant genomes have been uploaded in figshare (https://figshare.com/articles/dataset/220_genomes_input_file/29497892). Identification result files for these species are also in figshare (https://figshare.com/articles/dataset/CTDGs_result_of_220_plant_genomes/29497916). 

**For each genome, we prepared three input files for ctdgFinder:**  
**chromosomes.csv:** species,chromosome,length,bandwidth  
**genes_parsed.csv:** acc,species,chromosome,start,end,strand  
**seqs_hmmer.fa:** protein sequence of fasta format  

**Result of CTDGs for each species:**  
**Merged_CTDGs_genes_result_in_220_genomes.tar.gz** and **Merged_CTDGs_numbers_result_in_220_genomes.tar.gz** contain the final results of CTDGs identification in 220 plant genomes.  
***_merged_genes_result**: Gene ID, Species code, Chromosome ID, CTDG ID, member order  
***_merged_numbers_result**: Species code, Chromosome ID, CTDG ID, Member count  

**· Extract ctdgFinder results files: combines overlapping cluster results for all species listed in the files under CTDG_out directory and merges clusters with shared bases:**


```bash
$ python extract_result.py list ## get initial extracted results

## merge the "genes_result" of Aty, if not been merged completely, try 2-3 times until no "should have been merged" info showing
## round 1
$ python parse_result_of_CTDG_finder.py Aty
$ python validate_merge_result.py Aty
##round 2
$ mv Aty_processed_genes_results Aty_genes_result
$ python parse_result_of_CTDG_finder.py Aty
$ python validate_merge_result.py Aty
## round x
## ...

## we can also check out the partial overlapped information of the final merged genes_result
$ python overlap_analysis.py Aty
## finally we can generate the numbers_result file from the merged genes_file
$ python process_numbers_result.py Aty_processed_raw_genes_result Aty_processed_raw_numbers_result

```

**· Clustering rates analysis of different genes across 220 species**

Calculate the distribution matrix of clustering rates of genes of 7 traits in 220 species and obtain the protein sequences of all query
```bash
$ python analyze_plant_genes.py --result_dir ./ --clustered_dir ./ --species_order 220_spe_list_low2high_latin_abbr --gene_info 220_genes_parsed.csv --protein_file 220.pep --output_csv out_put/gene_distribution_stats.csv --output_seq_dir out_put
```

Select genes with the largest cluster change rate

```bash
$ Rscript ../custom_script/cluster_ratio_analyse.R -i ./gene_distribution_stats.csv -o ./cluster_ratio_analysis
```

Calculate the distribution matrix of clustering rates of different genes in different species and the bubble diagram of the first three traits

```bash
$ Rscript ../../custom_script/Bubble_plot_of_cluster_ratio.R -i gene_distribution_stats.csv -g group2spe.csv -n gene_family_name2_gene_ID.csv -o result
```





**· For each stress condition, caculate the DE-CTDGs according to the Mean Expression Difference and the Percentage of DEGs in Members, finally get the DE-CTDGs matrix for further analysis:**

```bash
## trimming
$ java -jar PATH/TO/trimmomatic-0.39.jar PE -threads 10 *.sra.fastq.gz *.sra.fastq.gz -baseout 1.clean_reads/* ILLUMINACLIP:PATH/TO/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

## quality control
fastqc -o . *.sra.fastq.gz

## align
$ hisat2 -t  -p 5 -x Aty -1 *.sra_1.fastq.gz -2 *.sra_2.fastq.gz -S *.sam
$ samtools view -S *.sam -b > *.bam
$ samtools sort -n *.bam -o *_nsorted.bam

## get reads count
$ htseq-count -n 10 -t gene -i ID -r name -f bam *_nsorted.bam Aty.gff3 > *.count

## calculate DE-CTDGs under each stress conditions for all size groups 
$ python identify_DECTDG.py --gene_file Aty_processed_raw_genes_result --number_file Aty_processed_raw_numbers_result --counts_file cold_counts_matrix.tsv --anno_file cold_pasAnno.csv --output_file cold_DEG_clusters.txt

## merge the DE-CTDGs lists under all conditions into the final matrix

$ python ~/01.project/2.207_ctdg/10.trancriptome_data_from_public/DE-CTDGs-v2.py --input_dir ./ --gff_file Aty.gff3 --gene_file Aty_processed_raw_genes_result --number_file Aty_processed_raw_numbers_result --output_dir DE-CTDG_output_no_sample --threads 10
# The input_dir must contain the listing files: gff3, _processed_raw_genes_result, _processed_raw_numbers_result, _counts_matrix.tsv (for each conditions), _pasAnno.csv (for each conditions)
```
**Now we have the DE-CTDGs bianry matrix, detailed results for each conditions and the TPM matrix. Let's move to analyze the functional enrichment for CTDGs in different size groups and DE-CTDGs of each single conditions**

  **For CTDGs in different size groups:**
  ```bash
  $ Rscript custom_script/GO_enrichment_for_diff_size_CTDGs.r
  ```

  **For DE-CTDGs of each single conditions:**
  ```bash
  $ Rscript ../../custom_script/GO_enrichment_for_multi_stress_DE-CTDGs-V2.r Aty
  ```


