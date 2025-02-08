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
**· For each stress condition, caculate the DE-CTDGs according to the Mean Expression Difference and the Percentage of DEGs in Members, finally get the DE-CTDGs matrix for further analysis:**
```bash
## calculate each conditon one by one and get the DE-CTDGs lists
python identify_DECTDG.py --gene_file Aty_processed_raw_genes_result --number_file Aty_processed_raw_numbers_result --counts_file cold_counts_matrix.tsv --anno_file cold_pasAnno.csv --output_file cold_DEG_clusters.txt
## merge the DE-CTDGs lists under all conditions into the final matrix
python create_expression_matrix.py all_CTDGs_list_of_Aty Aty_processed_raw_numbers_result expression_matrix.csv cold_DEG_clusters.txt drought_DEG_clusters.txt heat_DEG_clusters.txt light_DEG_clusters.txt salt_DEG_clusters.txt
```
