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
        # Rename chr[.] to [.] for matching with CADD
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
--regulatory \
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
# Note this needs to happen after combining with CADD annotations, otherwise bcftools annotation will only annotation the 1st match.
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
        # Fix a bug in VEP where gnomad AFs are assigned 'string' type. this line may not be necessary in newer VEP versions
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
# Note that for newer gencode version (> 37, possibly earlier versions) the 'lincRNA' gene_type annotation
# is removed and all lincRNA are categorized under lncRNA.
rule gencode:
    input:
        "data/gencode/{prefix}.gtf"
    output:
        gtf="data/gencode/{prefix}.chr-renamed.gtf",
        gene_len="data/gencode/{prefix}.chr-renamed.gene_lengths.tsv",
        bed="data/gencode/{prefix}.chr-renamed.bed",
        filt_bed="data/gencode/{prefix}.chr-renamed.protein_coding.lincRNA.bed",
        filt_bed_nov="data/gencode/{prefix}.chr-renamed.protein_coding.lincRNA.rm_ensemble_version.bed",
    conda: "envs/watershed.yml"
    shell:
        '''
        # Rename chromosome
        cat {input} | sed 's/^chr//g' > {output.gtf}
        # Calculate gene lengths for tpm
        gtftools -l {output.gene_len} {output.gtf}
        # Generate gene position bedfile 
        gtftools -g {output.bed} {output.gtf}
        # Filter by lincRNA / protein-coding convert to list of gene names
        cat {output.bed} | awk '{{if($7 == "lincRNA" || $7 == "protein_coding") {{print $0}} }}'  > {output.filt_bed}
        # Remove ENSEMBL version number (.*) for matching with VEP gene annotations.
        cat {output.filt_bed} | sed 's/\.[0-9]*//g' > {output.filt_bed_nov}
        '''

# Aggregate each individual's rare variants over each gene.
rule aggregate:
    input: 
        tsv="data/watershed/{prefix}.all.tsv",
        pc_linc=expand("data/gencode/{gencode}.protein_coding.lincRNA.rm_ensemble_version.bed", gencode=config["gencode"])
    output: "data/watershed/{prefix}.agg.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/agg.R {input.tsv} {input.pc_linc} {config[aggregate]} {config[types]} {config[vep_window]} > {output}
        '''   

# Add eOutliers and apply eOutlier filtering criteria.
rule add_eoutliers:
    input: 
        tsv="data/watershed/{prefix}.agg.tsv",
        scores="data/eoutliers/{prefix}.{cfg}.eOutliers.tsv"
    output:
        "data/watershed/{prefix}.{cfg}.eOutliers.agg.tsv"
    shell:
        '''
        Rscript scripts/add_eout.R {input.tsv} {input.scores} {config[pvalue]} {config[nout_std_thresh]} > {output}
        '''

# Add sOutliers and adjust minimum empirical pvalues.
rule add_soutliers:
    input: 
        tsv="data/watershed/{prefix}.{cfg}.agg.tsv",
        scores="data/soutliers/{prefix}.sOutliers.min_empirical_pvalue_gene.tsv"
    output:
        "data/watershed/{prefix}.{cfg}.sOutliers.agg.tsv"
    shell:
        '''
        Rscript scripts/add_sout.R {input.tsv} {input.scores} {config[pvalue]} {config[nout_std_thresh]} > {output}
        '''

# Label individuals with the same set of variants within each gene window.
rule label_pairs:
    input:
        "data/watershed/{prefix}.agg.tsv",
    output:
        "data/watershed/{prefix}.agg.pairlabel.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/label_pairs.R {input} > {output}
        ''' 

# Encode categorical variables (represented in a single column as comma-separated strings) as binary vectors.
rule encode_categorical:
    input: 
       tsv="data/watershed/{prefix}.agg.pairlabel.tsv",
       cat=config["categorical"]
    output: "data/watershed/{prefix}.agg.pairlabel.cat.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/encode_cat.R {input.tsv} {input.cat} > {output}
        '''

# Impute missing values
rule impute_missing:
    input:
       tsv="data/watershed/{prefix}.agg.pairlabel.cat.tsv",
       impute=config["impute"]
    output: "data/watershed/{prefix}.agg.pairlabel.cat.impute.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/impute_missing.R {input.tsv} {input.impute} {config[num_outliers]} > {output}
        '''

# Create test/train split
rule train_split:
    input: "data/watershed/{prefix}.agg.pairlabel.cat.impute.tsv",
    output: 
        tsv_test="data/watershed/{prefix}.agg.pairlabel.cat.impute.test.tsv",
        tsv_train="data/watershed/{prefix}.agg.pairlabel.cat.impute.train.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/train_split.R {input}
        '''

# Run Watershed
rule watershed:
    input: 
        full="data/watershed/{prefix}.agg.pairlabel.cat.impute.tsv",
        train="data/watershed/{prefix}.agg.pairlabel.cat.impute.train.tsv",
        test="data/watershed/{prefix}.agg.pairlabel.cat.impute.test.tsv",
    output: 
        eval="results/{prefix}/{seed}_evaluation_object.rds",
        predict="results/{prefix}/{seed}_posterior_probability.txt",
        cfg="results/{prefix}/{seed}_config.txt"
    conda: "envs/watershed.yml"
    shell:
        '''
        mkdir -p results/{wildcards.prefix}
        Rscript scripts/watershed.R {input.full} {input.train} {input.test} {config[num_outliers]} {config[pvalue]} {config[seed]} {config[C]} results/{wildcards.prefix}
        echo "p\t{config[pvalue]}" > {output.cfg}
        echo "C\t{config[C]}" >> {output.cfg}
        echo "nout_std_thresh\t{config[nout_std_thresh]}" >> {output.cfg}
        '''




