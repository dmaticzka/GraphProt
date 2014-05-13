################################################################################
### classification

## this is how I created the GraphProt test data for classification
#echo; echo; echo; echo prepare input data
#head -n 20 data/testclip.train.positives.fa > testclip.train.positives.fa
#head -n 20 data/testclip.train.negatives.fa > testclip.train.negatives.fa

echo; echo; echo; echo classification, cv
../GraphProt.pl -mode classification -action cv -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_cv --keep-tmp

echo; echo; echo; echo classification, train
../GraphProt.pl -mode classification -action train -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_train --keep-tmp

echo; echo; echo; echo classification, predict
../GraphProt.pl -mode classification -action predict -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_predict --keep-tmp

echo; echo; echo; echo classification, predict using only positives
../GraphProt.pl -mode classification -action predict    -fasta testclip.train.positives.fa -model CL_train.model  -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_predict_onlypos --keep-tmp

echo; echo; echo; echo classification nt margins
../GraphProt.pl -mode classification -action predict_profile --onlyseq -fasta testclip.train.positives.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_ntmargins --keep-tmp

echo; echo; echo; echo classification nt margins
../GraphProt.pl -mode classification -action predict_has --onlyseq -fasta testclip.train.positives.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -percentile 55 -prefix CL_has --keep-tmp

echo; echo; echo; echo classification motif
../GraphProt.pl -mode classification -action motif -fasta testclip.train.positives.fa -model CL_train.model  -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_motif --keep-tmp

echo; echo; echo; echo classification motif onlyseq
../GraphProt.pl -mode classification -action motif --onlyseq -fasta testclip.train.positives.fa -model CL_train.model -R 1 -D 0 -bitsize 10 -prefix CL_onlyseq --keep-tmp

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

################################################################################
### classification test for EDeN accuracy bug

echo; echo; echo; echo "train model for test-on-train results"

../GraphProt.pl --action train --fasta testedenacc.train.positives.fa --negfasta testedenacc.train.negatives.fa -prefix test_edenacc_train --keep-tmp

echo; echo; echo; echo "manual test-on-train"

../GraphProt.pl --action predict --fasta testedenacc.train.positives.fa --negfasta testedenacc.train.negatives.fa -model test_edenacc_train.model -prefix test_edenacc_manualtestontrain --keep-tmp