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
**· Extract ctdgFinder results files: combines overlapping cluster results for all species listed in the files under CTDG_out directory and merges clusters with shared bases:**

