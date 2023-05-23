library(WatershedR)



train="/oak/stanford/groups/smontgom/erobb/data/watershed/AF.all.ESN.hg38a.ID.ba.VEP.rare.agg.ws_train.tsv"
test="/oak/stanford/groups/smontgom/erobb/data/watershed/AF.all.ESN.hg38a.ID.ba.VEP.rare.agg.ws_test.tsv"
input = paste0("https://raw.githubusercontent.com/BennyStrobes/Watershed/",
     "master/example_data/watershed_example_data.txt")

# Run using Watershed approximate inference
evaluate_watershed(input_file = input,
                   model_name = "Watershed_approximate", 
                   number_of_dimensions = 3,
                   output_prefix = "watershed_approximate_n3")

