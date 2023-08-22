# watershed
Data preprocessing pipeline for [WatershedR](https://github.com/nicolerg/WatershedR).

## TODOs

High priority:
- [ ] Convert categorical variables to binary vector
- [ ] Add eOutlier scores
- [ ] Add pair labels

Mid priority:
- [ ] Gencode
    * [ ] Document where to download gencode files
    * [ ] Move content of scripts/gtf_filter.sh,pad.gtf.exons.py,process_gencode.sh into Snakemake pipeline
- [ ] Document configuration and all needed data input
- [ ] Add eOutlier residual calculation script to snakefile

Low priority:
- [ ] Add phyloP annotations
- [ ] Add distTSS, distTES annotations
- [ ] Add sOutlier scores

