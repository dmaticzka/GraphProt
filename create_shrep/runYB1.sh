perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/YB1_data_full_A.fa -t 1 -c 15 | gzip > YB1_data_full_A.shrep_t1_c15.gspan.gz
perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/YB1_data_full_B.fa -t 1 -c 15 | gzip > YB1_data_full_B.shrep_t1_c15.gspan.gz
perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/YB1_data_bruijn_A.fa -t 1 -c 15 | gzip > YB1_data_bruijn_A.shrep_t1_c15.gspan.gz
perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/YB1_data_bruijn_B.fa -t 1 -c 15 | gzip > YB1_data_bruijn_B.shrep_t1_c15.gspan.gz
