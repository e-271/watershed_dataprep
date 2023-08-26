
# Make symlinks to data folder and VEP folder from project root as follows:
"""
cd watershed
ln -s $DATA_PATH data
ln -s $VEP_PATH ensembl-vep
ln -s $VEP_HIDDEN .vep
"""


# TODO make a "rule all" and generalize other rule naming schemes
# TODO add all non-VCF inputs to config (CADD, LOFTEE path, Gencode, paths to configuration files)

configfile: "config/config.yaml"


# Annotate samples with CADD
# TODO may not be handling deletions properly due to differences in format between vcf and cadd
rule cadd:
    input:
        vcf="data/vcf/{prefix}.filt.ref_af.rare.vcf.gz",
        cadd_cols="config/cadd_columns",
        cadd="data/cadd/whole_genome_SNVs_inclAnno.tsv.gz",
        cadd_indel="data/cadd/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz"
    output: "data/vcf/{prefix}.filt.ref_af.rare.CADD.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
       '''
        sh scripts/anno_cadd.sh {input.cadd} {input.cadd_indel} {input.cadd_cols} {input.vcf} | bgzip > {output}
       tabix {output}
       '''

# Annotate rare variants with VEP
# Note that VEP does not support multithreading with the loftee plugin, so this needs to be run single-threaded.
# TODO VEP --distance=5000 by default for annotating upstream/downstream effects, we want this to be 10k presumably?
rule vep:
    input:
        vcf="data/vcf/{prefix}.filt.ref_af.rare.vcf.gz",
        loftee_path=".vep/Plugins"
    output: 
        vep="data/vcf/{prefix}.filt.ref_af.rare.VEP.vcf",
        vepgz="data/vcf/{prefix}.filt.ref_af.rare.VEP.vcf.gz",
    conda: "envs/vep.yml"
    shell:
        '''
        ensembl-vep/vep \
--verbose \
--vcf \
-i {input.vcf} \
-o {output.vep} \
--no_stats \
--force_overwrite \
--offline \
--dir_cache data/vep \
--plugin LoF,\
human_ancestor_fa:data/vep/hg38/human_ancestor.fa.gz,\
loftee_path:{input.loftee_path},\
conservation_file:data/vep/hg38/loftee.sql,\
gerp_bigwig:data/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
        bgzip --keep {output.vep}
        tabix {output.vepgz}
        '''

# Split all VEP annotations into INFO/* fields
rule split_vep:
    input: "data/vcf/{prefix}.filt.ref_af.rare.VEP.vcf.gz"
    output: "data/vcf/{prefix}.filt.ref_af.rare.VEP_split.vcf.gz"
    conda: "envs/watershed.yml"
    shell:
        '''
        fields=$(bcftools +split-vep {input} -l | cut -f 2 | tr '\n' ',')
        for l in {{1..22}} X Y; do echo chr$l $l; done > rename_chr.tmp
        bcftools +split-vep {input} -c $fields -d | bcftools annotate --rename-chrs rename_chr.tmp -x INFO/CSQ -o {output} 
        tabix {output}
        '''

# Can run VEP and CADD simultaneously & combine them afterwards to save time.
rule combine_annotations:
    input:
        vcf1="data/vcf/{prefix}.filt.ref_af.rare.CADD.vcf.gz",
        vcf2="data/vcf/{prefix}.filt.ref_af.rare.VEP_split.vcf.gz",
    output:
        "data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.vcf.gz",
    conda: "envs/watershed.yml"
    shell:
        '''
        bcftools annotate -a {input.vcf1} -c INFO {input.vcf2} -o {output}
        tabix {output}
        '''

# Split multi-sample VCF into a directory of single-sample VCFs
rule split_samples:
    input:
        "data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.vcf.gz"
    output:
        directory("data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.id_split")
    conda: "envs/watershed.yml"
    shell:
        '''
        bcftools +split {input} -o {output} -i 'GT="alt"'
        for f in {output}/*; do
            bgzip $f 
            tabix $f.gz
        done
        '''

# Filter & reformat gencode GTF file.
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
        # Filter & convert to (name, chr, start, end) tsv
        cat {output.genes} | awk '{{if($7 == "lincRNA" || $7 == "protein_coding") {{print $5 FS $1 FS $2 FS $3}} }}' > {output.tsv}
        '''

# Aggregate each individual's rare variants over each gene.
rule aggregate:
    input:
        vcf=directory("data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.id_split"),
        gencode=expand("data/gencode/{gencode}.gene_pos.protein_coding.lincRNA.tsv", gencode=config["gencode"]), 
        aggregate=config["aggregate"]
    output: 
        "data/watershed/{prefix}.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        sh scripts/agg.sh {input.vcf} {input.aggregate} {input.gencode} > {output}
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
        sh scripts/flatten.sh {input} {wildcards.type} > {output}
        '''

# Outlier scores are stored in an array of (gene) x (sample id).
rule add_outlier_scores:
    input:
        tsv="data/watershed/{prefix}.tsv",
        scores="data/outliers/{prefix}_{type}.split.tsv"
    output:
        tsv_tmp=temp("data/watershed/{prefix}.{type}_tmp.tsv"),
        tsv="data/watershed/{prefix}.{type}.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        # Join input file with outliers file
        xsv join Sample,Gene {input.tsv} Sample,Gene {input.scores} -o {output.tsv_tmp}
        # Remove duplicate columns used for matching
        xsv select '!Gene[1],Sample[1]' {output.tsv_tmp} -o {output.tsv} 
        '''   

# Label individuals with the same set of variants within each gene window.
rule label_pairs:
    input:
        "data/watershed/{prefix}.tsv",
    output:
        pairs=temp("data/watershed/{prefix}.pairs.tsv"),
        tsv_tmp=temp("data/watershed/{prefix}.pairlabel_tmp.tsv"),
        tsv="data/watershed/{prefix}.pairlabel.tsv"
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
       tsv="data/watershed/{prefix}.pairlabel.tsv",
       cat="config/categorical"
    output: "data/watershed/{prefix}.pairlabel.cat.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/encode_cat.R {input.tsv} {input.cat} > {output}
        '''

# Take norm of z-scores to make them comparable to p-values
rule normalize_zscores:
    input: "data/watershed/{prefix}.pairlabel.cat.tsv",
    output: "data/watershed/{prefix}.pairlabel.cat.normz.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/norm_zscores.R {input} {config[zscores]} > {output}
        '''   

# Impute missing values
rule impute_missing:
    input:
       vcf="data/watershed/{prefix}.pairlabel.cat.normz.tsv",
       impute="config/impute"
    output: "data/watershed/{prefix}.pairlabel.cat.normz.impute.tsv"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/impute_missing.R {input.vcf} {input.impute} > {output}
        '''

# Run Watershed
rule watershed:
    input: "data/watershed/{prefix}.pairlabel.cat.normz.impute.tsv",
    output: "data/watershed/{prefix}.pairlabel.cat.normz.impute_results"
    conda: "envs/watershed.yml"
    shell:
        '''
        Rscript scripts/watershed.R {input} {config[num_outliers]} > {output}
        '''




