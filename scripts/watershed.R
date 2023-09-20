library(WatershedR)

args = commandArgs(trailingOnly=TRUE)
input = args[1]
train = args[2]
test = args[3]
num_outliers = args[4]
out_dir = args[5]
out_prefix = tools::file_path_sans_ext(basename(input))
output=file.path(out_dir, out_prefix)

if (num_outliers==1) {
model_name = "RIVER"
} else {
# TODO when to use exact? depends on # dimenions I think?
model_name = "Watershed_approximate"
}

ws = evaluate_watershed(input_file = input,
                   model_name = model_name,
                   number_of_dimensions = as.numeric(num_outliers),
                   output_prefix = output,
                   n2_pair_pvalue_fraction = 0.01,
                   binary_pvalue_threshold = 0.01
           )

ws = predict_watershed(train, input,
                   number_dimensions = as.numeric(num_outliers),
                   model_name = model_name,
                   output_prefix = output,
                   binary_pvalue_threshold = 0.01
           )  





