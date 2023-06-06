devtools::install("/home/erobb/watershed/WatershedR")
library(WatershedR)

input="/oak/stanford/groups/smontgom/erobb/data/watershed/AF.all.ESN.hg38a.ID.ba.VEP.rare.agg.ws.pairs.tsv"
# Run using Watershed approximate inference
ws = evaluate_watershed(input_file = input,
                   model_name = "RIVER", 
                   number_of_dimensions = 1,
                   output_prefix = "watershed_approximate_n1")



print("Watershed AUC")
print(ws$auc[[1]]$evaROC$watershed_pr_auc)

print("GAM AUC")
print(ws$auc[[1]]$evaROC$GAM_pr_auc)


names(ws$auc[[1]]$evaROC)

# list of length 18
names(ws$model_params)

# Model parameters, 
# edges between outliers (dim=1 for us)
"theta_pair"           
"theta_singleton"      

# feat_dim x outlier_dim
"theta"               

# Means per gene/ind, very large
"mu"                   
"mu_pairwise"       

# eOutlier posterior, length $num_samples
"posterior"           
"posterior_pairwise"   

# Feature metadata
"num_samples"          
"num_genomic_features"
"number_of_dimensions" 

# ??
"phi"                  

# Lambda hyperparameter, selected by k-fold cross-validation
# This is regularization weight
"lambda"              
"lambda_singleton"     
"lambda_pair"          

# ??
"pseudoc"             

# VI hyperparams
"vi_step_size"         
"vi_thresh"            

"model_name"  




# example dim=3
#input="https://raw.githubusercontent.com/BennyStrobes/Watershed/master/example_data/watershed_example_data.txt"

# example dim=1
#input="https://raw.githubusercontent.com/BennyStrobes/Watershed/v1.0.0/example_data/river_example_data_pheno_1.txt"

