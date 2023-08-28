# watershed
Data preprocessing pipeline for [WatershedR](https://github.com/nicolerg/WatershedR).

## Set up CADD
Download CADD files to match your genome build, and place in the folder `data/cadd`. 

CADD (HG38):
- SNVs
  - [tsv.gz](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz) (313GB)
  - [tbi](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz.tbi)
- Indels
  - [tsv.gz](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz) (7.6GB)
  - [tbi](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz.tbi)

If you use a different genome build or location for CADD, update the  `cadd` and `cadd_indel` fields in `config/config.yaml`. You can also configure the annotation columns and aggregation functions for CADD by editing the `config/cadd_columns` file.

## Set up VEP

1. Create a VEP conda environment. Make sure the match the VEP version to your Ensembl version:

    `conda env create --name vep -c bioconda ensembl-vep=$ENSEMBL_VERSION`

2. Download cache files for VEP annotation and place in the folder `data/vep`. Run `tar xzf` to unzip this after downloading:
    - [Ensembl 110 / GRCh38](https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz) (13GB)
    - [Ensembl 110 / GRCh37](https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/#:~:text=homo_sapiens_vep_110_GRCh37.tar.gz) (20GB)
    - Other versions can to be located within [https://ftp.ensembl.org/pub](https://ftp.ensembl.org/pub) or can be installed using VEP's `install.PL` script.

For full VEP installation instructions see the [VEP documentation](http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html).

## Set up Loftee

1. Install the [Loftee](https://github.com/konradjk/loftee) plugin for VEP in the folder `.vep/Plugins`:
    - GRCh38: `git clone https://github.com/konradjk/loftee --branch grch38 .vep/Plugins`
    - GRCh37: `git clone https://github.com/konradjk/loftee .vep/Plugins`

2. Download Loftee files to `data/vep`:
    - GERP (GRCh38):
      - [bw](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw)
    - Human ancestor (GRCh38):
      - [fa.gz](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz)
      - [fai](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai)
      - [gzi](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi)
    - PhyloCSV (GRCh38):
      - [sql.gz](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz) (unzip after downloading)

## Run the pipeline

Example inputs can be found in `data/vcf/EXAMPLE.vcf.gz` and `data/outliers/EXAMPLE_eOutliers.tsv`.
1. Install Watershed conda environment:

    `conda install --file envs/watershed.yml`

2. Run the pipeline (including CADD & VEP annotation, data preparation, & Watershed training):

    `snakemake --cores $N data/watershed/EXAMPLE.eOutliers.pairlabel.cat.normz.impute.format_results`

The full SnakeMake workflow is shown below:

![DAG](docs/dag.svg?raw=true)

## Notes on configuration
- `config/config.yaml` defines filepaths for annotation files and other general settings for the Snakemake pipeline.
- `config/cadd_columns` defines which CADD annotations to use, as well as the aggregation function for variants returning multiple CADD annotations (e.g. those with multiple variant categories). Available aggregation functions are `max`, `unique`, `min`, `append`.
- `config/categorical` defines VEP and CADD annotations which need to be converted to categorical variables.
- `config/aggregate` defines the aggregation functions over each gene region. Available aggregation functions are `max`, `unique`, `min`, `append`.
- `config/impute` defines the imputation values used to fill any missing annotations.


## Notes on speed

- The `VEP`, `CADD`, and `aggregate` steps can take 24+ hours on a dataset with ~10M rare variants. 
- The `aggregate` step can be highly parallelized.
- Watershed training only takes a few minutes.
- VEP with the Loftee plugin does not support multi-threading, as their SQL queries are not thread-safe.

# TODOs

High priority:
- [ ] Add example data files
    * [ ] VCF input (with only minimal required fields, like INFO/AF and quality)
    * [ ] Outlier scores file (to match rule add_outlier_scores input.scores)
- [ ] Delete extra individuals for N>2 pairs
- [ ] Update `config/categorical` and `scripts/encode_cat.R` to define 'keep' categories instead of 'drop' categories
- [ ] Match all configuration files to original Watershed paper: [Ferraro & Strober et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7646251/) supplementary table S3.

Mid priority:
- [ ] Add phyloP annotations
- [ ] Add distTSS, distTES annotations
- [ ] Add MAF and num_rare_variants to annotations

