
# Make symlinks to data folder and VEP folder from project root as follows:
"""
cd watershed
ln -s $DATA_PATH data
ln -s $VEP_PATH ensembl-vep
ln -s $VEP_HIDDEN .vep
"""


# TODO add conda envs
# TODO make a "rule all" and generalize other rule naming schemes
# TODO add all non-VCF inputs to config
# TODO put config 

configfile: "config/config.yaml"

# Filter to VQSR "PASS"
rule filter:
    input: "data/vcf/{prefix}.vcf.gz"
    output: "data/vcf/{prefix}.filt.vcf.gz"
    shell:
       '''
       bcftools view -f PASS {input} -o {output}
       tabix {output}
       '''

# Update AF field in VCF using AC / AN
rule update_afs:
    input: "data/vcf/{prefix}.filt.vcf.gz"
    output: "data/vcf/{prefix}.filt.sample_af.vcf.gz"
    shell:
        '''
        bcftools +fill-tags {input} -o {output} -- -t AF
        tabix {output}
        '''

# Update population AF using a reference AF
rule update_pop_afs:
    input:
       vcf="data/vcf/{prefix}.filt.vcf.gz",
       ref=expand("data/vcf/{reference}.vcf.gz", reference=config["reference"])
    output: "data/vcf/{prefix}.filt.ref_af.vcf.gz"
    shell:
        '''
        bcftools annotate -a {input.ref} -c INFO/AF {input.vcf} -o {output}
        tabix {output}
        '''

# Filter positions by MAF<0.01, AC>0. Also splits multi-allelic positions into biallelic records.
rule filter_rare:
    input: "data/vcf/{prefix}.filt.ref_af.vcf.gz" 
    output: "data/vcf/{prefix}.filt.ref_af.rare.vcf.gz" 
    shell:
        '''
        bcftools norm -m -any {input} |  bcftools view --exclude "AF>0.01 | AC=0" -o {output}
        tabix {output}
        '''

# Annotate samples with CADD
# TODO may not be handling deletions properly due to differences in format between vcf and cadd
rule cadd:
    input:
        vcf="data/vcf/{prefix}.filt.ref_af.rare.vcf.gz",
        cadd_cols="data/cadd/cadd_columns",
        cadd="data/cadd/whole_genome_SNVs_inclAnno.tsv.gz",
        cadd_indel="data/cadd/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz"
    output: "data/vcf/{prefix}.filt.ref_af.rare.CADD.vcf.gz"
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
    shell:
        '''
        bcftools annotate -a {input.vcf1} -c INFO {input.vcf2} -o {output}
        tabix {output}
        '''

rule split_samples:
    input:
        "data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.vcf.gz"
    output:
        directory("data/vcf/{prefix}.filt.ref_af.rare.CADD.VEP_split.id_split") # Is a directory of <id>.vcf
    shell:
        '''
        bcftools +split {input} -o {output} -i 'GT="alt"'
        for f in {output}/*; do
            bgzip $f 
            tabix $f.gz
        done
        '''

# input.gencode contains 4 columns: gene_id, chr, pos_start, pos_end
# input.fields contains vcf INFO fields to include in Watershed annotations (one per line)
# TODO add gencode file to config 
# Very rough estimate is that this will take about 24 hours (guessing 5s / gene). 
rule aggregate:
    input:
        vcf=directory("data/vcf/{prefix}"),
        gencode="data/gencode/gencode.v43.gene_pos.protein_lincRNA.tsv",
        fields="config/watershed_fields"
    output: 
        "data/watershed/{prefix}.agg.tsv"
    shell:
        '''
        sh scripts/agg.sh {input.vcf} {input.fields} {input.gencode} > {output}
        '''

# I am not sure this is something bcftools can do.
# It may be possible but is not elegant. I think it's best if I can append together all categoricals during the aggregation step, and then handle max/min/convert to categorical after it is in a TSV format.
rule encode_categorical:
    input: 
       vcf="data/vcf/{prefix}.vcf",
       cat="config/categories" # TODO maybe this comes from the config file
    output: "data/vcf/{prefix}.cat.vcf"
    shell:
        '''
        '''


