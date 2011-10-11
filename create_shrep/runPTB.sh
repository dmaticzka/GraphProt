# steffens parameters
perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_A.fa -t 3 -c 20 -M 5 | gzip > PTB_data_full_A.shrep_t3_c20_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_B.fa -t 3 -c 20 -M 5 | gzip > PTB_data_full_B.shrep_t3_c20_M5.gspan.gz

# probability cutoff using shape types
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_A.fa -t 1 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_A.shrep_t1_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_B.fa -t 1 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_B.shrep_t1_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_A.fa -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_A.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_B.fa -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_B.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_A.fa -t 5 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_A.shrep_t5_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/full/PTB_data_full_B.fa -t 5 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_full_B.shrep_t5_T005_M5.gspan.gz

# probability cutoff with PTB CLIP-seq
perl RNAtools/fasta2shrep_gspan_T.pl -fa PTB_pos_pad25.target100.fa -t 3 -M 5 -wins 50 -shift 10 | gzip > PTB_pos_pad25.target100.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa PTB_pos_pad25.target100.fa -t 3 -M 5 -wins 50 -shift 10 -cue | gzip > PTB_pos_pad25.target100.shrep_t3_T005_M5_cue.gspan.gz

# probability cutoff using shape types and weakly structured sequences
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_A.fa -t 1 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_A.shrep_t1_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_B.fa -t 1 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_B.shrep_t1_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_A.fa -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_A.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_B.fa -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_B.shrep_t3_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_A.fa -t 5 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_A.shrep_t5_T005_M5.gspan.gz
perl RNAtools/fasta2shrep_gspan_T.pl -fa ../data/RNAcompete_derived_datasets/fastas/weak/PTB_data_bruijn_B.fa -t 5 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_data_bruijn_B.shrep_t5_T005_M5.gspan.gz

# compute shreps for final PTB experiment
perl RNAtools/fasta2shrep_gspan.pl -fa PTB_train.fasta -wins "50,100" -shift 25 -cue -c 20 -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > PTB_train.gspan.gz
perl RNAtools/fasta2shrep_gspan.pl -fa ANXA7.fasta -wins "50,100" -shift 25 -cue -c 20 -t 3 -M 5 | awk -f removeEmptyLines.awk | gzip > ANXA7.fasta.gz
