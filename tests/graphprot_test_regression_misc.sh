################################################################################
### regression

## this is how I created the GraphProt test data for regression
#echo; echo; echo; echo prepare input data
#cp data/test_data_full_A.train.fa ./
#grep '^>' test_data_full_A.train.fa | awk '{print $2}' > test_data_full_A.train.affys

export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH

echo; echo; echo; echo regression, cv
export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
../GraphProt.pl -mode regression -action cv -fasta test_data_full_A.train.fa -affinities test_data_full_A.train.affys -prefix REG_cv -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp

echo; echo; echo; echo regression, train
export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
../GraphProt.pl -mode regression -action train -fasta test_data_full_A.train.fa  -affinities test_data_full_A.train.affys -prefix REG_train -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp

echo; echo; echo; echo regression, predict
export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
../GraphProt.pl -mode regression -action predict -fasta test_data_full_A.train.fa -affinities test_data_full_A.train.affys -model REG_train.model -prefix REG_predict -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp

echo; echo; echo; echo regression, predict no affys
export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
../GraphProt.pl -mode regression -action predict -fasta test_data_full_A.train.fa -model REG_train.model -prefix REG_predict_noaffys -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp

