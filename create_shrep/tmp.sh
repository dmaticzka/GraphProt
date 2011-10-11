perl RNAtools/fasta2shrep_gspan_T.pl -fa PTB_pos_pad25.target100.fa -t 3 -M 5 -wins 50 -shift 10 | gzip > PTB_pos_pad25.target100.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa PTB_pos_pad25.target100.fa -t 3 -M 5 -wins 50 -shift 10 -cue | gzip > PTB_pos_pad25.target100.shrep_t3_T005_M5_cue.gspan.gz
