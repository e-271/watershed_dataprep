#devtools::install("/home/erobb/WatershedR")
library(WatershedR)

args = commandArgs(trailingOnly=TRUE)
input = args[1]
num_outliers = args[2]
output = args[3]

if (num_outliers==1) {
model_name = "RIVER"
} else {
# TODO when to use exact? depends on # dimenions I think?
model_name = "Watershed_approximate"
}


ws = evaluate_watershed(input_file = input,
                   model_name = model_name,
                   number_of_dimensions = num_outliers,
                   output_prefix = output
                   )







