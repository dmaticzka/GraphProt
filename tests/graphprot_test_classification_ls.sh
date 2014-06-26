echo; echo; echo; echo classification, ls
../GraphProt.pl -mode classification -action ls -fasta testclip.train.positives.fa -negfasta testclip.train.negatives.fa -prefix CL_ls --keep-tmp

echo; echo; echo; echo train model using parameters file from linesearch
../GraphProt.pl -mode classification -action train -fasta testclip.train.positives.fa -negfasta testclip.train.negatives.fa -params CL_ls.params -prefix CL_train_from_ls --keep-tmp
