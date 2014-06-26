echo; echo; echo; echo regression, ls
export PATH=/home/maticzkd/src/libsvm-3.12/:$PATH
../GraphProt.pl -mode regression -action ls -fasta test_data_full_A.train.fa -affinities test_data_full_A.train.affys -prefix REG_ls -abstraction 1 -R 0 -D 0 -epsilon 0.11 -c 11 -bitsize 10 --keep-tmp

echo; echo; echo; echo train model using parameters file from linesearch
../GraphProt.pl -mode regression -action train -fasta test_data_full_A.train.fa -affinities test_data_full_A.train.affys -params REG_ls.params -prefix REG_train_from_ls --keep-tmp
