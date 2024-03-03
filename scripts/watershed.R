#devtools::install("/home/erobb/WatershedR")
library(WatershedR)

args = commandArgs(trailingOnly=TRUE)
input = args[1]
train = args[2]
test = args[3]
num_outliers = args[4]
pvalue = as.numeric(args[5])
seed = as.numeric(args[6])
C = as.numeric(args[7])
out_dir = args[8]

out_prefix = toString(seed)
output=file.path(out_dir, out_prefix)

# compute true outlier proportion for n2_pair_pvalue_fraction
df = read.table(input,header=T)
pvalue_prop = sum(abs(df$eOutliers) < pvalue) / dim(df)[1]
print(sprintf("total measured pvalue proportion %.4f", pvalue_prop) )

if (num_outliers==1) {
model_name = "RIVER"
} else {
model_name = "Watershed_exact"
}

set.seed(seed)
ws = evaluate_watershed(input_file = input,
                   model_name = model_name,
                   number_of_dimensions = as.numeric(num_outliers),
                   output_prefix = output,
                   n2_pair_pvalue_fraction = pvalue_prop,
                   binary_pvalue_threshold = pvalue,
                   dirichlet_prior = C,
                   l2_prior_parameter=0.01 # note: default is set too high in WatershedR package
           )

ws = predict_watershed(train, input,
                   number_dimensions = as.numeric(num_outliers),
                   model_name = model_name,
                   output_prefix = output,
                   binary_pvalue_threshold = pvalue,
                   dirichlet_prior_parameter = C,
                   l2_prior_parameter=0.01
           )  





