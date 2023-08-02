

pop=$1
#if [ $pop == "YRI" ] 
#then
#af=0.025
#else
#af=0.01
#fi

echo $pop #$af


#python scripts/eqtls.sh $1
#python preprocess/residuals.R $1
python preprocess/make_dataset_mp.py --pop $1
python preprocess/add_gencode.py --pop $1
python preprocess/add_phylop.py --pop $1
python preprocess/aggregate_genes.py --pop $1
python preprocess/add_eout.py --pop $1
#python preprocess/add_sout.py --pop $1
python preprocess/add_pair_labels.py --pop $1 

