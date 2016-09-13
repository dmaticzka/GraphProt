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

echo; echo; echo; echo classification nt margins high affinity sites
../GraphProt.pl -mode classification -action predict_has --onlyseq -fasta testclip.train.positives.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -percentile 55 -prefix CL_has --keep-tmp

echo; echo; echo; echo classification motif
../GraphProt.pl -mode classification -action motif -fasta testclip.train.positives.fa -model CL_train.model  -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix CL_motif --keep-tmp
gunzip -f CL_motif.motif.gspan.gz

echo; echo; echo; echo classification motif onlyseq
../GraphProt.pl -mode classification -action motif --onlyseq -fasta testclip.train.positives.fa -model CL_train.model -R 1 -D 0 -bitsize 10 -prefix CL_onlyseq --keep-tmp

### test profile calculation for structure models and sequences exactly 300nt
## all viewpoints

echo; echo; echo; echo classification nt margins allcaps
../GraphProt.pl -mode classification -action predict_profile -fasta testclip.twoseqs300nt_allcaps.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix testclip.twoseqs300nt_allcaps.margins --keep-tmp

## only center viewpoints
echo; echo; echo; echo classification nt margins
../GraphProt.pl -mode classification -action predict_profile -fasta testclip.twoseqs300nt.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix testclip.twoseqs300nt.margins --keep-tmp

## only center viewpoints motif
echo; echo; echo; echo classification nt margins
../GraphProt.pl -mode classification -action motif -fasta testclip.twoseqs300nt.fa -model CL_train.model -abstraction 1 -R 0 -D 0 -epochs 2 -lambda 0.1 -bitsize 10 -prefix testclip.twoseqs300nt.motif --keep-tmp

gunzip -f *.gz
