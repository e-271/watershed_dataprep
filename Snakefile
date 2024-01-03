""" Prepare Watershed input data from a vcf. """
configfile: "config/config.yaml"

# Filter positions by MAF<0.01, AC>0, and <10% missing sample genotypes
# Splits multi-allelic record into biallelic records (e.g. if ALT=T,A split into 2 records with ALT=T, ALT=A).
rule filter_rare:
    input: "data/vcf/{prefix}.vcf.gz"
    output: "data/vcf/{prefix}.snp.rare.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Rename chr[.] to [.]
        for l in {{1..22}} X Y; do echo chr$l $l; done > rename_chr.tmp
        bcftools norm -m -any {input} | 
            bcftools view  -f PASS -i 'TYPE=\"snp\" & AF<=0.01 & AC>0 & F_MISSING<0.1' |
            bcftools annotate --rename-chrs rename_chr.tmp -o {output} 
        tabix {output}
        '''

# Annotate samples with CADD
rule cadd:
    input: "data/vcf/{prefix}.rare.vcf.gz",
    output: 
        vcf="data/vcf/{prefix}.rare.CADD.vcf",
        gz="data/vcf/{prefix}.rare.CADD.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
       '''
       sh scripts/anno_cadd.sh {config[cadd]} {config[cadd_cols]} {input} > {output.vcf}
       bgzip --keep {output.vcf}
       tabix {output.gz}
       '''

# Annotate rare variants with VEP
# VEP does not support multithreading with the loftee plugin, so this needs to be run single-threaded.
# The --custom flag format varies by VEP version. For version-specific documentation see the VEP archives: 
# http://useast.ensembl.org/info/website/archives/index.html
# For gnomad v4 use the jointly called AFs: AF_joint_afr, AF_joint_amr, ...
# VEP90:
#--custom {config[gnomad]},gnomADg,vcf,exact,0,AF_joint_afr,AF_joint_amr,AF_joint_asj,AF_joint_eas,AF_joint_sas,AF_joint_fin,AF_joint_nfe \
rule vep:
    input: "data/vcf/{prefix}.rare.vcf.gz"
    output: 
        vep="data/vcf/{prefix}.rare.VEP.gnomad.vcf",
        vepgz="data/vcf/{prefix}.rare.VEP.gnomad.vcf.gz",
    conda: "envs/watershed.yml"
    shell:
        '''
        vep \
--verbose \
--vcf \
-i {input} \
-o {output.vep} \
--distance {config[vep_window]} \
--no_stats \
--force_overwrite \
--offline \
--dir_cache data/vep \
--dir_plugins {config[vep_plugins_dir]} \
--custom file={config[gnomad]},short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe \
--plugin LoF,\
human_ancestor_fa:{config[human_ancestor]},\
loftee_path:{config[vep_plugins_dir]},\
conservation_file:{config[conservation_file]},\
gerp_bigwig:{config[gerp_bigwig]}
        bgzip --keep {output.vep}
        tabix {output.vepgz}
    '''

# Combine CADD and VEP annotated vcfs.
rule combine_annotations:
    input:
        vcf1="data/vcf/{prefix}.rare.CADD.vcf.gz",
        vcf2="data/vcf/{prefix}.rare.VEP.gnomad.vcf.gz",
    output:
        "data/vcf/{prefix}.rare.CADD.VEP.gnomad.vcf.gz",
    conda: "envs/watershed.yml"
    shell:
        '''
        bcftools annotate -a {input.vcf1} -c INFO {input.vcf2} -o {output}
        tabix {output}
        '''

# Split all VEP annotations into INFO/* fields, and split multiple gene annotations into separate records
rule split_vep:
    input: "data/vcf/{prefix}.rare.CADD.VEP.gnomad.vcf.gz"
    output: "data/vcf/{prefix}.rare.CADD.VEP.gnomad.split.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        fields=$(bcftools +split-vep {input} -l | cut -f 2 | tr '\n' ',')
        bcftools +split-vep {input} -c $fields -d | \
            bcftools annotate -x INFO/CSQ -o {output}
        tabix {output}
        '''

# Filter to variants that are rare in all gnomAD populations
rule rare_gnomad:
    input: "data/vcf/{prefix}.rare.CADD.VEP.gnomad.split.vcf.gz"
    output: 
      rh="data/vcf/{prefix}.rare.CADD.VEP.gnomad_rh.split.vcf.gz",
      rare="data/vcf/{prefix}.rare.CADD.VEP.gnomad_rh_rare.split.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        # fixes a bug in VEP where gnomad AFs are assigned 'string' type. this line may not be necessary in newer VEP versions
        bcftools head {input} | sed '/^##INFO=<ID=gnomADg_AF/s/Type=String/Type=Float/g' | bcftools reheader -h - {input} -o {output.rh}
        tabix {output.rh}
        bcftools view -i 'gnomADg_AF_joint_afr<=0.01 & gnomADg_AF_joint_amr<=0.01 & gnomADg_AF_joint_asj<=0.01 & gnomADg_AF_joint_eas<=0.01 & gnomADg_AF_joint_fin<=0.01 & gnomADg_AF_joint_nfe<=0.01 & gnomADg_AF_joint_sas<=0.01' {output.rh} -o {output.rare}
        tabix {output}
        '''

# Format rare variants to a tsv
rule tsv_format:
    input: "data/vcf/{prefix}.rare.CADD.VEP.gnomad_rh_rare.split.vcf.gz"
    output: "data/watershed/{prefix}.all.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        sh scripts/format.sh {input} {config[format]} > {output}
        '''

# Format gencode file to use in gene-level aggregation. Filters by protein-coding / lincRNA genes.
rule gencode:
    input:
        expand("data/gencode/{prefix}.gtf", prefix=config["gencode"])
    output:
        bed="data/gencode/{prefix}.bed",
        filt_bed="data/gencode/{prefix}.protein_coding.lincRNA.bed",
        tsv="data/gencode/{prefix}.protein_coding.lincRNA.gene_names.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Generate gene position bedfile 
        gtftools -g {output.bed} {input}
        # Filter, convert to list of gene names, and remove version number (.*)
        cat {output.bed} | awk '{{if($7 == "lincRNA" || $7 == "protein_coding") {{print $0}} }}' | sed 's/\.[0-9]*//g' > {output.filt_bed}
        # List of gene names only
        cut -f5 {output.filt_bed} > {output.tsv}
        '''

# Split outlier scores to 1 line per sample, and remove .[0-9]* in Ensemble identifiers
rule format_outliers:
    input:
        "data/outliers/{prefix}_{type}-{cfg}.tsv"
    output:
        "data/outliers/{prefix}_{type}-{cfg}.split.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Flatten expression matrix to (gene,subject) rows
        sh scripts/flatten.sh {input} {wildcards.type} > {output}
        '''

# Aggregate each individual's rare variants over each gene.
rule aggregate:
    input: 
        tsv="data/watershed/{prefix}.all.tsv",
        pc_linc=expand("data/gencode/{gencode}.protein_coding.lincRNA.bed", gencode=config["gencode"])
    output: "data/watershed/{prefix}.agg.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/agg.R {input.tsv} {input.pc_linc} {config[aggregate]} {config[types]} {config[vep_window]} > {output}
        '''

# Add outlier scores to tsv. 
rule add_outlier_scores:
    input:
        tsv="data/watershed/{prefix}.agg.tsv",
        scores="data/outliers/{prefix}_{type}.split.tsv"
    output:
        tsv_tmp=temp("data/watershed/{prefix}.{type}_tmp.tsv"),
        tsv="data/watershed/{prefix}.{type}.agg.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Join input file with outliers file
        xsv join Sample,Gene {input.tsv} Sample,Gene {input.scores} -o {output.tsv_tmp}
        # Remove duplicate columns used for matching
        xsv select '!Gene[1],Sample[1]' {output.tsv_tmp} -o {output.tsv} 
        '''   

# Filter to protein-coding / lincRNA genes with at least 1 eOutlier
rule filter_genes:
    input: 
        tsv="data/watershed/{prefix}.eOutliers{cfg}.agg.tsv",
    output:
        "data/watershed/{prefix}.eOutliers{cfg}.agg.filt.tsv"
    shell:
        '''
        # Filter tsv to these genes
        Rscript scripts/filter_genes.R {input.tsv} {config[zthreshold]} {config[nout_std_thresh]} > {output}
        '''

# Label individuals with the same set of variants within each gene window.
rule label_pairs:
    input:
        "data/watershed/{prefix}.agg.filt.tsv",
    output:
        "data/watershed/{prefix}.agg.filt.pairlabel.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/label_pairs.R {input} > {output}
        ''' 

# Encode categorical variables (represented in a single column as comma-separated strings) as binary vectors.
rule encode_categorical:
    input: 
       tsv="data/watershed/{prefix}.agg.filt.pairlabel.tsv",
       cat=config["categorical"]
    output: "data/watershed/{prefix}.agg.filt.pairlabel.cat.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/encode_cat.R {input.tsv} {input.cat} > {output}
        '''

# Impute missing values
rule impute_missing:
    input:
       tsv="data/watershed/{prefix}.agg.filt.pairlabel.cat.tsv",
       impute=config["impute"]
    output: "data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        uscript scripts/impute_missing.R {input.tsv} {input.impute} {config[num_outliers]} > {output}
        '''

# Drop and rename columns for Watershed
rule format:
    input:
       tsv="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.tsv",
       rename=config["rename"]
    output: 
        tsv="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.tsv",
        tsv_test="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.test.tsv",
        tsv_train="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.train.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/train_split.R {input.tsv} {input.rename}
        '''

# Run Watershed
rule watershed:
    input: 
        full="data/watershed/{prefix}.pairlabel.cat.impute.format.tsv",
        train="data/watershed/{prefix}.pairlabel.cat.impute.format.train.tsv",
        test="data/watershed/{prefix}.pairlabel.cat.impute.format.test.tsv",
    output: 
        eval="results/watershed/{seed}/{prefix}.pairlabel.cat.impute.format_evaluation_object.rds",
        predict="results/watershed/{seed}/{prefix}.pairlabel.cat.impute.format_posterior_probability.txt"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/watershed.R {input.full} {input.train} {input.test} {config[num_outliers]} {config[pvalue]} {config[zthreshold]} {wildcards.seed} results/watershed/{wildcards.seed}
        '''




