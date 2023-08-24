#devtools::install("/home/erobb/WatershedR")
library(WatershedR)

args = commandArgs(trailingOnly=TRUE)
pop = args[1]
s = as.numeric(args[2])
data_dir = "/oak/stanford/groups/smontgom/erobb/data"


print("starting...")
print(pop)
print(s)
set.seed(s)
input = sprintf("%s/watershed/all.%s.30x.ID.VEP.bedtools.rare.ws.gencode.phyloP.agg.eout.pairs.tsv", data_dir, pop)
#input = sprintf("%s/watershed/AF.all.%s.hg38a.ID.ba.VEP.gencode.phyloP.agg.eout.pairs.rare.ws.tsv", data_dir, pop)
#input = sprintf("/Volumes/T7/watershed/AF.all.%s.hg38a.ID.ba.VEP.gencode.phyloP-241.agg.eout.pairs.rare.ws.tsv", pop)
ws = evaluate_watershed(input_file = input,
                   model_name = "RIVER",
                   number_of_dimensions = 1,
                   output_prefix = sprintf("RIVER_n1_30x_%s_seed%s", pop, s),
                    #output_prefix = sprintf("RIVER_n1_%s_phyloP-241_seed%s", pop, s),
                   )







