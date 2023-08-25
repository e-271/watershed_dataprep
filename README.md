# watershed
Data preprocessing pipeline for [WatershedR](https://github.com/nicolerg/WatershedR).

## TODOs

High priority:
- [x] Convert z-scores to p-scores
- [x] Impute null values
- [ ] Prepare test set (delete extra individuals for N>2 pairs)
- [ ] Update script to run Watershed 
- [x] Convert categorical variables to binary vector
- [x] Add eOutlier scores
- [x] Add pair labels

Mid priority:
- [ ] Gencode
    * [ ] Document where to download gencode files
    * [ ] Move content of scripts/gtf_filter.sh,pad.gtf.exons.py,process_gencode.sh into Snakemake pipeline
- [ ] Document configuration and all needed data input
- [ ] Add eOutlier residual calculation script to snakefile
- [ ] Review/standardize Snakemake naming scheme
- [ ] Add example data files
    * [ ] VCF input (with only minimal required fields, like INFO/AF and quality)
    * [ ] Outlier scores file (to match rule add_outlier_scores input.scores)
    * [ ] Watershed input

Low priority:
- [ ] Add phyloP annotations
- [ ] Add distTSS, distTES annotations
- [ ] Add sOutlier scores

