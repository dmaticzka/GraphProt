# classification, ls
GraphProt_distribution/GraphProt.pl -mode classification -action ls    -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa

# classification, cv
GraphProt_distribution/GraphProt.pl -mode classification -action cv    -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa -R 1 -D 0 -bitsize 10

# classification, train
GraphProt_distribution/GraphProt.pl -mode classification -action train -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa

# classification, predict
GraphProt_distribution/GraphProt.pl -mode classification -action predict    -fasta testclip.train.positives.fa  -negfasta testclip.train.negatives.fa -model GraphProt.model -R 1 -D 0 -bitsize 10

# classification, predict using only positives
GraphProt_distribution/GraphProt.pl -mode classification -action predict    -fasta testclip.train.positives.fa -model GraphProt.model -R 1 -D 0 -bitsize 10

# classification nt margins
GraphProt_distribution/GraphProt.pl -mode classification -action predict_profile --onlyseq -fasta testclip.train.positives.fa -model GraphProt.model -R 1 -D 0 -bitsize 10

# classification motif
GraphProt_distribution/GraphProt.pl -mode classification -action motif -fasta testclip.train.positives.fa -model GraphProt.model -R 1 -D 0 -bitsize 10

# classification motif onlyseq
GraphProt_distribution/GraphProt.pl -mode classification -action motif --onlyseq -fasta testclip.train.positives.fa -model GraphProt.model -R 1 -D 0 -bitsize 10

# regression, ls
GraphProt_distribution/GraphProt.pl -mode regression -action ls    -fasta test_data_full_A.test.fa  -affinities test_data_full_A.test.affys

# regression, cv
GraphProt_distribution/GraphProt.pl -mode regression -action cv    -fasta test_data_full_A.test.fa  -affinities test_data_full_A.test.affys -R 1 -D 0 -bitsize 10

# regression, train
GraphProt_distribution/GraphProt.pl -mode regression -action train -fasta test_data_full_A.test.fa  -affinities test_data_full_A.test.affys

# regression, predict
GraphProt_distribution/GraphProt.pl -mode regression -action predict    -fasta test_data_full_A.test.fa  -affinities test_data_full_A.test.affys -model GraphProt.model -R 1 -D 0 -bitsize 10
