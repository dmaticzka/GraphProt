################################################################################
### classification test for EDeN accuracy bug

echo; echo; echo; echo "train model for test-on-train results"

../GraphProt.pl --action train --fasta testedenacc.train.positives.fa --negfasta testedenacc.train.negatives.fa -prefix test_edenacc_train --keep-tmp

echo; echo; echo; echo "manual test-on-train"

../GraphProt.pl --action predict --fasta testedenacc.train.positives.fa --negfasta testedenacc.train.negatives.fa -model test_edenacc_train.model -prefix test_edenacc_manualtestontrain --keep-tmp
