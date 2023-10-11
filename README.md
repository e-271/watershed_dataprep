# watershed
Data preprocessing pipeline for [WatershedR](https://github.com/nicolerg/WatershedR).

## Download CADD files
1. Download CADD files to match your genome build, and place in the folder `data/cadd`.

    - SNVs (HG38):
      - [tsv.gz](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz) (313G)
      - [tbi](https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz.tbi)

If you use a different genome build or filepath for CADD, update `config/config.yaml`.

## Set up VEP

1. Create a VEP conda environment. Make sure the match the VEP version to your Ensembl version:

    `conda env create --name vep -c bioconda ensembl-vep=$ENSEMBL_VERSION`

2. Download cache files for VEP annotation and place in the folder `data/vep`. Run `tar xzf` to unzip this after downloading:
    - [Ensembl 110 / GRCh38](https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz) (13G)
    - [Ensembl 110 / GRCh37](https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/#:~:text=homo_sapiens_vep_110_GRCh37.tar.gz) (20G)
    - Other versions can to be located within [https://ftp.ensembl.org/pub](https://ftp.ensembl.org/pub) or can be installed using VEP's `install.PL` script.

3. Set the `${VEP_PLUGINS_DIR}` path variable. By default VEP installs plugins in `$HOME/.vep/Plugins`, but Conda environments will store this in `$CONDA_PREFIX/share/ensembl-vep-$VERSION`.
    - If you use `ensemble-vep` conda packages from the `dnachun` channel, the `${VEP_PLUGINS_DIR}` environment variable will be set to the correct location by default.

For full VEP installation instructions see the [VEP documentation](http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html).

If you use a different genome build, Ensembl version or filepath for VEP, update `config/config.yaml`.

## Set up Loftee

VEP installs a version of loftee for HG37 among the default plugins. To use it you will need to install cache files (step 2). If your input is aligned to HG38 you can use one of the `ensemble-vep` conda packages from the `dnachun` channel which have loftee-hg38 installed, or copy the loftee-hg38 perl files from Github (step 1).


1. [If using HG38] Install the [Loftee (hg38)](https://github.com/konradjk/loftee) plugin for VEP:
    - GRCh38: `git clone https://github.com/konradjk/loftee --branch grch38 $VEP_PLUGINS_DIR`

2. Download Loftee files to `data/vep`:
    - GERP (GRCh38):
      - [bw](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw) (12G)
    - Human ancestor (GRCh38):
      - [fa.gz](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz) (844M)
      - [fai](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai)
      - [gzi](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi)
    - PhyloCSV:
      - [sql.gz](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz) (29M) unzip after downloading

If you use a different genome build or filepath for Loftee, update `config/config.yaml`.

## Download Gencode file

1. Download the GTF file corresponding to your build, unzip it, and and place in `data/gencode`:
    - [Gencode v44 / GRCh38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) (47M)
    - [Gencode v43 / GRCh37](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/gencode.v43lift37.annotation.gtf.gz) (62M)

Update the gencode version in `config/config.yaml`.

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

