
# Make symlinks to data folder and VEP folder from project root as follows:
"""
cd watershed
ln -s $DATA_PATH data
ln -s $VEP_PATH ensembl-vep
ln -s $VEP_HIDDEN .vep
"""


# TODO make a "rule all" and generalize other rule naming schemes

configfile: "config/config.yaml"

# Filter positions by MAF<0.01, AC>0. 
# Also splits multi-allelic record into biallelic records (e.g. if ALT=T,A split into 2 records with ALT=T, ALT=A).
rule filter_rare:
    input: "data/vcf/{prefix}.vcf.gz"
    output: "data/vcf/{prefix}.snp.rare.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Rename chr[.] to [.]
        for l in {{1..22}} X Y; do echo chr$l $l; done > rename_chr.tmp
        bcftools norm -m -any {input} | 
            bcftools view -i 'TYPE=\"snp\" & AF<=0.01 & AC>0' | 
            bcftools annotate --rename-chrs rename_chr.tmp -o {output} 
        tabix {output}
        '''

# Annotate samples with CADD
# TODO remove cadd_indels
rule cadd:
    input: "data/vcf/{prefix}.rare.vcf.gz",
    output: "data/vcf/{prefix}.rare.CADD.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
       '''
        sh scripts/anno_cadd.sh {config[cadd]} {config[cadd_indels]} {config[cadd_cols]} {input} | bgzip > {output}
       tabix {output}
       '''

# Annotate rare variants with VEP
# Note that VEP does not support multithreading with the loftee plugin, so this needs to be run single-threaded.
rule vep:
    input: "data/vcf/{prefix}.rare.vcf.gz"
    output: 
        vep="data/vcf/{prefix}.rare.VEP.vcf",
        vepgz="data/vcf/{prefix}.rare.VEP.vcf.gz",
    conda: "envs/vep.yml"
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
--plugin LoF,\
human_ancestor_fa:{config[human_ancestor]},\
loftee_path:{config[loftee_path]},\
conservation_file:{config[conservation_file]},\
gerp_bigwig:{config[gerp_bigwig]}
        bgzip --keep {output.vep}
        tabix {output.vepgz}
    '''

# Can run VEP and CADD simultaneously & combine them afterwards to save time.
# This needs to happen before vep-split due to bcftools annotate not annotating multiple matching lines.
rule combine_annotations:
    input:
        vcf1="data/vcf/{prefix}.rare.CADD.vcf.gz",
        vcf2="data/vcf/{prefix}.rare.VEP.vcf.gz",
    output:
        "data/vcf/{prefix}.rare.CADD.VEP.vcf.gz",
    conda: "envs/watershed.yml"
    shell:
        '''
        bcftools annotate -a {input.vcf1} -c INFO {input.vcf2} -o {output}
        tabix {output}
        '''

# Split all VEP annotations into INFO/* fields, and split multiple gene annotations into separate records
rule split_vep:
    input: "data/vcf/{prefix}.rare.CADD.VEP.vcf.gz"
    output: "data/vcf/{prefix}.rare.CADD.VEP_split.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        fields=$(bcftools +split-vep {input} -l | cut -f 2 | tr '\n' ',')
        bcftools +split-vep {input} -c $fields -d | \
            bcftools annotate -x INFO/CSQ -o {output}
        tabix {output}
        '''

# Format rare variants to a tsv
rule tsv_format:
    input: "data/vcf/{prefix}.rare.CADD.VEP_split.vcf.gz"
    output: "data/watershed/{prefix}.all.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        sh scripts/format.sh {input} {config[aggregate]} > {output}
        '''

# Format gencode GTF file & filter by protein-coding / lincRNA genes.
rule gencode:
    input:
        expand("data/gencode/{prefix}.gtf", prefix=config["gencode"])
    output:
        genes="data/gencode/{prefix}.genes.bed",
        tsv="data/gencode/{prefix}.gene_pos.protein_coding.lincRNA.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Generate gene position bedfile 
        gtftools -g {output.genes} {input}
        # Filter, convert to list of gene names, and remove version number (.*)
        cat {output.genes} | awk '{{if($7 == "lincRNA" || $7 == "protein_coding") {{print $5}} }}' | sed 's/\.[0-9]*$//g' > {output.tsv}
        '''

# Split outlier scores to 1 line per sample, and remove .[0-9]* in Ensemble identifiers
rule format_outliers:
    input:
        "data/outliers/{prefix}_{type}.tsv"
    output:
        "data/outliers/{prefix}_{type}.split.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Flatten expression matrix to (gene,subject) rows
        sh scripts/flatten.sh {input} {wildcards.type} > {output}
        '''

# Aggregate each individual's rare variants over each gene.
# The R version requires a lot of memory but is much faster. To reduce memory, split the inputs by subject and/or chromosome.
# There is also a bcftools version in agg.sh, which takes several days to run but uses little memory.
rule aggregate:
    input: "data/watershed/{prefix}.all.tsv"
    output: "data/watershed/{prefix}.agg.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/agg.R {input} {config[aggregate]} {config[types]} > {output}
        '''

# Outlier scores are stored in an array of (gene) x (sample id).
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
        tsv="data/watershed/{prefix}.eOutliers.agg.tsv",
        pc_linc=expand("data/gencode/{gencode}.gene_pos.protein_coding.lincRNA.tsv", gencode=config["gencode"]) 
    output:
        "data/watershed/{prefix}.eOutliers.agg.filt.tsv"
    shell:
        '''
        # Filter tsv to these genes
        Rscript scripts/filter_genes.R {input.tsv} {input.pc_linc} {config[zthreshold]} {config[maxnout]} > {output}
        '''

# Label individuals with the same set of variants within each gene window.
rule label_pairs:
    input:
        "data/watershed/{prefix}.agg.filt.tsv",
    output:
        pairs=temp("data/watershed/{prefix}.agg.filt.pairs.tsv"),
        tsv_tmp=temp("data/watershed/{prefix}.agg.filt.pairlabel_tmp.tsv"),
        tsv="data/watershed/{prefix}.agg.filt.pairlabel.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Get unique (gene,position,alt) combinations and label each
        echo -e "Gene\\tPOS\\tALT\\tPair\\tPairN" > {output.pairs}
        tail -n +2 {input} | \
            cut -f 2-4 | sort | uniq -c -d | nl | \
            awk '{{print $3 "\\t" $4 "\\t" $5 "\\t" $1 "\\t" $2}}' \
            >> {output.pairs}

        # Join pair labels with input file
        xsv join --left Gene,POS,ALT {input} Gene,POS,ALT {output.pairs} -o {output.tsv_tmp}
        # Remove duplicate columns used for matching, and delete POS & ALT fields 
        xsv select '!Gene[1],POS[0],POS[1],ALT[0],ALT[1]' {output.tsv_tmp} -o {output.tsv} 
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
        Rscript scripts/impute_missing.R {input.tsv} {input.impute} {config[num_outliers]} > {output}
        '''

# Drop and rename columns for Watershed
rule format:
    input:
       tsv="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.tsv",
       drop=config["drop"],
       rename=config["rename"]
    output: 
        tsv_tmp=temp("data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format_tmp.tsv"),
        tsv="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.tsv",
        tsv_test="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.test.tsv",
        tsv_train="data/watershed/{prefix}.agg.filt.pairlabel.cat.impute.format.train.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # drop columns
        drop='!'
        while read drop_col; do drop="${{drop}}${{drop_col}}[0],"; done < {input.drop}
        xsv select $drop {input.tsv} -o {output.tsv_tmp}
 
        # rename columns
        head -n 1 {output.tsv_tmp} | \
            awk '{{if(NR==FNR){{ a[$1]=$2 }} else{{ for(i=1;i<=NF;i++) {{ $i=a[$i]?a[$i]:$i }} print $0 }}}}' \
            FS=',' {input.rename} FS='\\t' - > {output.tsv}
        tail -n +2 {output.tsv_tmp} >> {output.tsv}
     
        # create 'train' and 'test' sets for predict_watershed function
        NF=$(head -n 1 {output.tsv}  | sed 's/\s/\\n/g' | wc -l)
        head -n 1 {output.tsv} > {output.tsv_test}
        cat {output.tsv} |  awk -v nf=$NF '$nf != "NA"' >> {output.tsv_test}
        head -n 1 {output.tsv} > {output.tsv_train}
        cat {output.tsv} |  awk -v nf=$NF '$nf == "NA"' >> {output.tsv_train}
        '''

# Run Watershed
rule watershed:
    input: 
        full="data/watershed/{prefix}.pairlabel.cat.impute.format.tsv",
        train="data/watershed/{prefix}.pairlabel.cat.impute.format.train.tsv",
        test="data/watershed/{prefix}.pairlabel.cat.impute.format.test.tsv",
    output: 
        eval="results/watershed/{prefix}.pairlabel.cat.impute.format_evaluation_object.rds",
        predict="results/watershed/{prefix}.pairlabel.cat.impute.format_posterior_probability.txt"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/watershed.R {input.full} {input.train} {input.test} {config[num_outliers]} results/watershed
        '''




