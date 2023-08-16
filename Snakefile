
# Make symlinks to data folder and VEP folder from project root as follows:
"""
cd watershed
ln -s $DATA_PATH data
ln -s $VEP_PATH ensembl-vep
ln -s $VEP_HIDDEN .vep
"""


# TODO add conda env

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
# Due to a bug in vcfanno, we need to combine the CADD SNV and indel annotation TSV files before running this.
rule cadd:
    input:
        vcf="data/vep/{prefix}.filt.ref_af.rare.vcf.gz",
        config="config/cadd.conf"
    output: "data/vep/{prefix}.filt.ref_af.rare.CADD.vcf",
    shell:
       '''
       bcftools view {input.vcf} | sed 's/chr//' | vcfanno {input.config} > {output}
       '''

# Annotate rare variants with VEP
# Note that VEP does not support multithreading with the loftee plugin, so this needs to be run single-threaded.
# TODO setup a Singularity image with loftee and put on github
rule vep:
    input:
        vcf="data/vcf/{prefix}.filt.ref_af.rare.vcf.gz",
        loftee_path=".vep/Plugins"
    output: "data/vep/{prefix}.filt.ref_af.rare.VEP.vcf"
    shell:
        '''
ensembl-vep/vep \
--verbose \
--vcf \
-i {input.vcf} \
-o {output} \
--no_stats \
--force_overwrite \
--offline \
--dir_cache data/vep \
--plugin LoF,\
human_ancestor_fa:data/vep/hg38/human_ancestor.fa.gz,\
loftee_path:{input.loftee_path},\
conservation_file:data/vep/hg38/loftee.sql,\
gerp_bigwig:data/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
        '''

# Convert to Watershed format
# Split VEP annotations and subset to 'consequence' and 'LoF' fields
# Filter samples by minimum AF<0.01 and print 1 sample per line
rule watershed_format:
    input: "data/vep/{prefix}.ref_af.rare.VEP.CADD.vcf"
    output: "data/vep/{prefix}.ref_af.rare.VEP.CADD.split.vcf"
    shell:
        '''
         bcftools +split-vep {input} -f '[%SAMPLE %CHROM %POS %REF %ALT %TYPE %GT %AF %Allele %Gene %Consequence %LoF\n]' -i 'GT="alt" & AF<0.01' > {output}
        '''


