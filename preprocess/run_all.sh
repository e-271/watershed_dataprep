

pop=$1

python make_dataset_mp.py --pop $1
python add_gencode.py --pop $1
python add_phylop.py --pop $1
python aggregate_dataset.py --pop $1
python add_pair_labels.py --pop $1 

