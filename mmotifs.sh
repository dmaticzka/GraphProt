# remake logos for all proteins with finished computations
for f in *.train.model
do
PREFIX=`basename $f .train.model`
	# make sure that sequence parts are recalculated properly
	echo rm -rf $PREFIX.test.*_top_wins*
	# make logos
	echo make \
	$PREFIX.test.sequence_top_wins.truncated.logo.png \
	$PREFIX.test.struct_annot_top_wins.truncated.logo.png \
	$PREFIX.test.struct_annot_top_wins.truncated.pup.logo.png
done
