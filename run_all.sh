for f in *.train.model
do
FPRE=`basename $f .train.model`
echo "echo $FPRE"
# calculate nt margins of all shreps
echo make \
$FPRE.test.sequence_top_wins.truncated.logo.png \
$FPRE.test.struct_annot_top_wins.truncated.logo.png \
-e MARGINS_WINDOW=12 -e MARGINS_MEASURE=mean -e TOP_WINDOWS=1000 -e GRAPH_TYPE=CONTEXTSHREP -e DO_LINESEARCH=YES -e LINESEARCH_INPUT_SIZE=1000 -e DO_SGDOPT=NO -o $FPRE.param
done
