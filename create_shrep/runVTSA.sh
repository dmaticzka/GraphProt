time perl RNAtools/fasta2shrep_gspan.pl -fa ../seqs/fastas/full/VTS1_data_full_A.fa -t 3 | gzip > VTS1_data_full_A.shrep_t3.gspan.gz
time perl RNAtools/fasta2shrep_gspan.pl -fa ../seqs/fastas/full/VTS1_data_full_B.fa -t 3 | gzip > VTS1_data_full_B.shrep_t3.gspan.gz

time perl RNAtools/fasta2shrep_gspan.pl -fa ../seqs/fastas/full/VTS1_data_full_A.fa -t 3 -e 10 > VTS1_data_full_A.shrep_t3_e10.gspan
time perl RNAtools/fasta2shrep_gspan.pl -fa ../seqs/fastas/full/VTS1_data_full_B.fa -t 3 -e 10 > VTS1_data_full_B.shrep_t3_e10.gspan

time perl RNAtools/fasta2shrep_gspan.pl -fa ../seqs/fastas/full/VTS1_data_full_A.fa -t 3 -e 20 > VTS1_data_full_A.shrep_t3_e20.gspan
