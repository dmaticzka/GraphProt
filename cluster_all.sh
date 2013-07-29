for f in *.train.model
do
FPRE=`basename $f .train.model`
echo "echo $FPRE"
# select some sequences
echo "cat $FPRE.train.positives.fa | /usr/local/perl/bin/perl /usr/local/user/RNAtools/fastapl -p1 | head -n 4000 > $FPRE.test.positives.fa"
echo "touch $FPRE.test.negatives.fa"
# calculate nt margins of all shreps
echo "echo 'make $FPRE.test.nt_margins $FPRE.test.struct_annot.nt_subset $FPRE.test.sequence_top_wins $FPRE.test.struct_annot_top_wins $FPRE.test.struct_annot_top_wins.pup -e MARGINS_WINDOW=12 -e MARGINS_MEASURE=mean -e TOP_WINDOWS=1000 -e GRAPH_TYPE=CONTEXTSHREP -e DO_LINESEARCH=YES -e LINESEARCH_INPUT_SIZE=1000 -e DO_SGDOPT=NO -o $FPRE.param' | qsub -m ea -M maticzkd@informatik.uni-freiburg.de -cwd -notify -e '\$JOB_NAME.jid\$JOB_ID.sgeout' -o '\$JOB_NAME.jid\$JOB_ID.sgeout' -l h_vmem=10G -N XXX_$FPRE -v OMP_NUM_THREADS=1 -v OMP_STACKSIZE=20M"
echo "sleep 2"
done
