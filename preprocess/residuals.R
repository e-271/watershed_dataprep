#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
pop = args[1]
#pops = c("ESN", "GWD", "LWK", "MKK", "YRI", "MSL")


data_dir="/oak/stanford/groups/smontgom/erobb/data"

files = c("ESN.SVs.11.5.txt",
        "GWD.SVs.13.6.txt",
        "LWK.SVs.10.4.txt",
        "MKK.SVs.14.4.txt",
        "YRI.SVs.8.2.txt",
        "MSL.SVs.11.4.txt")
names(files) <- c("ESN", "GWD", "LWK", "MKK", "YRI", "MSL")

# Read in deseq2 expression data
expr_file=sprintf('%s/deseq/deseq2ScaledGeneCounts.%s.bed',data_dir,pop)
expr = read.delim(expr_file)
row.names(expr) = expr[,"ID"] # set row name to gene ID

# Read in covariates matrix
covs_file = sprintf("/oak/stanford/groups/smontgom/mdegorte/durga/africa/fastqtl/SVs/%s/%s",pop,files[[pop]])
covs = read.table(covs_file, header = T, sep = '\t', row.names = 1)

# Read in top eQTL
eqtl_file=sprintf('%s/eqtls/sorted.byGene.dist.hwe.af.%s.top1.eQTL.nominal.hg38a.txt', data_dir, pop)
eqtl_df = read.table(eqtl_file, header = T, row.names = 'feature')
# Read in vcf data for each top eQTL
vcf_file = sprintf("%s/eqtls/AF.all.%s.hg38aID.eQTLs.ba.vcf", data_dir, pop)
vcf_df = read.table(vcf_file, header=T, skip=0, row.names="GENE")

# Filter to IDs with available covariates
id_keep = colnames(expr)[5:length(colnames(expr))]
id_keep = intersect(id_keep,colnames(covs))
id_keep = intersect(id_keep,colnames(vcf_df))
covs = covs[,id_keep]
expr = expr[,id_keep]

# Subset to protein-coding & lincRNA genes
pc = eqtl_df["geneType"] == "protein_coding"
linc = eqtl_df["geneType"] == "lincRNA"
eqtl_df = eqtl_df[pc | linc,]

# 1 if alt allele is present, 0 otherwise
vcf_df[,id_keep] = sapply(vcf_df[,id_keep] == '0|0', as.numeric)
vcf_df = vcf_df[,id_keep]

# Scale to mean 0 variance 1 (over each gene)
expr = scale(t(expr))

# # For each gene in the expression file, perform a linear regression 
resids = matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)

for(i in 1:ncol(expr)){ # for each gene! seems slow but not too bad
    gene = colnames(expr)[i]
    if (!all(is.na(vcf_df[gene,]))) {
        gcovs = rbind(covs, vcf_df[1,])
    } else {
        gcovs = covs
    }
    data = as.data.frame(cbind(expr[, i], t(gcovs)))
    colnames(data) = c('deseq2ScaledGeneCounts', rownames(gcovs))
    model = lm(deseq2ScaledGeneCounts ~ ., data = data)
    resids[, i] = model$residuals
}
# Center and scale, then transpose
resids = t(scale(resids))

# Write out the residuals
out_file = sprintf('%s/eoutliers/%s_exprResiduals.tsv', data_dir, pop)
write.table(matrix(c('Id', colnames(resids)), nrow = 1), out_file, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(resids, out_file, row.names = T, col.names = F, quote = F, sep = '\t', append = T)


#####
# Label outliers

QTHRESH = F

if (QTHRESH) {
    q = quantile(resids)
    z_thresh = as.numeric(abs(q["75%"] - q["25%"]))
} else {
    z_thresh = 3
}

# AFGR individuals have more on average, so we can choose a higher threshold than Watershed
z_ind_filt = 75
r_file = sprintf('%s/eoutliers/%s_exprResiduals.tsv', data_dir, pop)

# Write out the residuals
resids = read.table(r_file, header=T, row.names=1)
resids = data.frame(t(resids))

w = abs(resids) > z_thresh
w = data.frame(sapply(data.frame(w), function(x) as.integer(x)))
rownames(w) = rownames(resids)
indiv_noutliers = data.frame(rowSums(w))
indiv_keep = rownames(indiv_noutliers)[which(indiv_noutliers < z_ind_filt)]

browser()

mean(indiv_noutliers[,])
sqrt(var(indiv_noutliers[,]))
median(indiv_noutliers[,])
print(sprintf("Keeping %i / %i individuals with z_thresh = %#.2f, max outliers = %i",length(indiv_keep), length(rownames(w)), z_thresh, z_ind_filt))
write.table(indiv_keep, sprintf('%s/eoutliers/ids_outlier_filtered_%s_t%#.2ff%i.txt', data_dir, pop, z_thresh, z_ind_filt), quote = F, sep = '\t', col.names = F, row.names = F)
write.table(w, sprintf('%s/eoutliers/eOutlier_scores_%s_t%#.2f.txt', data_dir, pop, z_thresh))




