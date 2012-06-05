package StructureLibrary::Structure;
use Cwd 'abs_path';
use strict;
use StructureLibrary::RNAtoolsConfig; # Config.pm must be in the same directory

# global variables
# This declares the named variables as package globals in the current package.
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %CONFIG);
*CONFIG = \%StructureLibrary::RNAtoolsConfig::CONFIG; # Global Configuration now aliased to global %CONFIG


require Exporter;
#require Tools;
#require Sequence;
our @ISA = qw(Exporter);

#our $VERSION = '0.01';
#
## Package global variables
#our $DOTBRACKET = "().-";
## matching brackets ( ): base i pairs base j
## . : base is unpaired
## - : structure unkown
#
#our $CONSTRAINT = "().|x<>";
## | : paired with another base
## . : no constraint at all
## x : base must not pair
## < : base i is paired with a base j<i
## > : base i is paired with a base j>i
## matching brackets ( ): base i pairs base j
#
#our %STR_FORMATS;
#$STR_FORMATS{DOTBRACKET} = $DOTBRACKET;
#$STR_FORMATS{CONSTRAINT} = $CONSTRAINT;

our @EXPORT = qw(
			accuracy
			call_RNAfold_parse_output
			call_RNAfold
			check_bps_wrt_seq
			structure2basepairs
			check_str_format
			ensemble_energies_and_probability_with_constraint
			parse_dotplot_return_bp_probs
			parse_dotplot_return_bp_probs_and_sequence
			parse_PUfile_get1u
			parse_PUfile_getallu
			parse_PUfile_getsomeu
			prob_base_i_is_paired
			read_seq_structure_fasta_file
			convert_structure_format
			remove_non_compatible_bps_wrt_seq
			basepairs2dotbracket
			add_dotplot_to_fabrizio_graph
			convert_dotplot_to_fabrizio_graph
			call_rfold
			rfold_1seq_and_modify_output
			convert_bp_to_accfile
			read_rfold_prob_file_into_basepairs_hash
			);
our @EXPORT_OK = qw(
				CONFIG
				convert_dotplot_to_fabrizio_graph
				add_dotplot_to_fabrizio_graph
				
			);
	#			%CONFIG		
##################################################################################
##################################################################################
## Package Structure.pm 	AUTHOR(S) = Sita Lange
## This package includes a collection of methods that are useful when handling 
## the secondary structure of RNA molecules, especially in combination with the
##  Vienna Package.
##################################################################################
##################################################################################
			
			
##################################################################################
# Calculates the accuracy of the given structure within a reference sequence, 
# although the structure can also be across the entire length of the reference
# sequence. The accuracy is calculated according to the given dotplot, i.e. base-
# pair probabilities and the unpaired probabilities, which are either given or
# calculated from the base-pair probabilities. The unpaired probabilities are
# calculated by summing over all base-pair probabilities and taking 1 minus this
# sum. Of course this only works with real probabilities. Some of the local folding
# methods do not generate real probabilities, so this is then the wrong calculation.
# The accuracy score consists of the
# sum of all base-pair probabilities for each base-pair within the given structure
# multiplied by two and added to that is the sum of all unpaired probabilites for
# all unpaired bases in the given structure. If you are unsure of the structure
# of any bases these can be ignored. The notation of the structures must be given
# in either alphabet defined in %STR_FORMATS. The default is the dotbracket format
# in $DOTBRACKTET.
# accuracy = sum_over_i:j_in_struct(P(i,j) x 2) 
#				+ sum_over_i_unpaired_in_struct(PU(i_unpaired)) 
# NOTE: ALL INDICES COUNT FROM 1 TO N, except of course within arrays. This is due
#		to the counting in the Vienna package files.
#
# Input: 
#		struct		The structure string in the one of the %STR_FORMATS alphabets.
#		str_f		(Optional) Format of structure in struct; one of %STR_FORMATS. 
#					Default is $DOTBRACKET, the dot-bracket format.
#		start		This is the starting position of the structure in the sequence.
#					(counting 1 to n)
#		bp_prob		A hash of ALL base-pair probabilities for the reference sequence,
#					where the key is "i:j" and i < j and they correspond to the given 
#					sequence.
#		pu			(Optional) An array with ALL unpaired probabilities for the 
#					the reference sequence (same length as bp_prob)
#					If this option is not given, the unpaired probabilites
#					are calculated from the base-pair probabilites as 
#					PU(i) = 1-sum_over_all_j(P(i,j)). NOTE: this only works for 
#					the global probabilites of RNAfold and NOT for RNAplfold!
#					You can parse PUfiles using parse_PUfile_get1u($PUfile,1,$num_header)
#		pufile		(Optional) A PUfile can be given instead of pu, so that the
#					unpaired probabilities can be read directly from the file.
#		num_header	(Compulsory if pufile is given) The number of header lines in the pufile.
#		seq			(Optional) The original/reference sequence. 
#					The structure is checked to see if it feasible given
#					this sequence (i.e. can the base-pairs form and are the
#					indices within the correct range?). 
#		seq_len		(Compulsory) Give the reference sequence length 
#					Also the indices of
#					the PU values are checked to be within the sequence boundaries.
#		normalise	(Optional) Normalise the accuracy score according to the number of
#					probabilities in the sum
#		silent		(Optional) If you are very sure of what you are doing you can suppress all
#					outputs, such as warnings.
# Output:
#		An array containing
#		(1) The accuracy score of the complete structure
#		(2) The accuracy score of the base-pairs
#		(3) The accuracy score of the unpaired bases
#		If the normalise flag, then each of these scores are normalised, which
#		means divided by the number of bases that go into the accuracy calculation.
#
# Note this method has been checked against an implementation from Kousik/Steffen
# and should be correct.
##################################################################################
sub accuracy{
	my ($struct, $str_f, $start, $bp_prob_href, $pu_aref, $pufile, $num_header, 
		$seq, $seq_len, $normalise, $silent) = @_;
	my $FUNCTION = "accuracy in Structure.pm";
	
	# check input
	die ("INPUT ERROR in $FUNCTION:\nWhere does the structure start in you sequence?\n") unless ($start);
	die ("INPUT ERROR in $FUNCTION:\nWhere is the input structure?\n") unless ($struct);
	die ("INPUT ERROR in $FUNCTION:\nWhere are the base-pair probabilites?\n") unless ($bp_prob_href);
	print STDERR ("WARNING in $FUNCTION:\nNo unpaired probabilities given, are you using real probabilities?\n") unless ($pu_aref || $pufile||$silent);
	die ("INPUT ERROR in $FUNCTION:\nYou can not give both pu and pufile. Please try again!\n") if ($pu_aref && $pufile);
	if($seq){
		$seq_len = length($seq) unless ($seq_len);
	} else {
		die ("INPUT ERROR in $FUNCTION:\nThe reference sequence length ".
			"was not given!\n") unless ($seq_len);
	}
	
	## check that the structure can still fit into the sequence, wrt its start pos
	my $str_len = length($struct);
	my $str_end = $start + $str_len - 1;
	if($str_end > $seq_len){
		die ("INPUT ERROR in $FUNCTION: the structure cannot fit into the given".
			"sequence w.r.t the given starting position!\n");
	}
#	print "seq_len:$seq_len str len: $str_len\t$struct with start $start\n";
	
	# accuracy score
	my $acc = 0;
#	my $normalise_sum = 0;
	my $normalise_bpacc = 0;
	my $normalise_puacc = 0;
	my $bpacc = 0;
	my $puacc = 0;
	
	# get base-pairs from structure	counting 1 to n
	$str_f = $CONFIG{STR_FORMATS}->{DOTBRACKT} unless ($str_f);
	my $bp_href = structure2basepairs($struct, $start, $str_f);
	
	# check structure
	check_bps_wrt_seq($seq, $bp_href) if($seq);
	
	## cannot have no base-pairs, this wouldn't make sense
	my @bp_str = keys(%{$bp_href});
	die "ERROR in $FUNCTION: this is not a valid structure! $struct\n" 
		if(@bp_str < 2);
	
	# add the probability of each base-pair in the structure and times it by 2
	foreach my $bp (@bp_str){
		if(exists($bp_prob_href->{$bp})){
			$bpacc += $bp_prob_href->{$bp} * 2;
		}
		$normalise_bpacc += 2;
	}
	
	# get unpaired positions relative to the reference sequence, counting 0 to n-1
	my @struct_a = split("", $struct);
	my @unpaired = ();
	# get unpaired symbol
	my $unp_symb = "";
	$unp_symb = "." if($str_f eq $CONFIG{STR_FORMATS}->{DOTBRACKET});
	$unp_symb = "x" if($str_f eq $CONFIG{STR_FORMATS}->{CONSTRAINT});
	die "ERROR in $FUNCTION: The structure format is unkown, i.e. the unpaired symbol\n" 
		unless ($unp_symb);
	my $n = @struct_a;
	for (my $i = 1 ; $i <= $n ; $i ++ ){
		push (@unpaired, ($i+$start-1)) if ($struct_a[$i-1] eq $unp_symb);
	}
	
	# get unpaired probabilities if not given
	# unpaired probabilities given as a file
	if($pufile) {
#		($num_header) or die ("ERROR in $FUNCTION: The number of header lines has not been given!\n");
		my @pu = parse_PUfile_get1u($pufile, 1, $num_header);
		$pu_aref = \@pu;
				
	# no unpaired probabilities are given
	} else {
		my $paired_i_href = prob_base_i_is_paired($bp_prob_href);
		$pu_aref = ();
		
		for (my $i = 1 ; $i <= $seq_len ; $i ++ ){
			push(@{$pu_aref}, (1-$paired_i_href->{$i}));
		}
	}
	
	# final check
	my $pu_len = @{$pu_aref};
	if($pu_len != $seq_len){
		die ("ERROR in $FUNCTION:\n The PU values do not match the length of the reference sequence!\n");
	}
	
	# add unpaired probabilities
	foreach (@unpaired){
		$puacc += $pu_aref->[$_ - 1];
		++ $normalise_puacc;
	}
	
	die ("ERROR in $FUNCTION: No structure probabilities have been added to the accuracy score!\n") 
		unless ($normalise_bpacc && $normalise_puacc);
	
	$acc = $bpacc + $puacc;
	$acc = $acc/($normalise_bpacc+$normalise_puacc) if ($normalise);
	$bpacc = $bpacc/$normalise_bpacc if ($normalise);
	$puacc = $puacc/$normalise_puacc if ($normalise);
#	print STDERR "base pair contribution: $bpacc\nunpaired contribution: $puacc\nfinal acc: $acc\n";
	
#	print "TEST ACCURACY: 1:$acc 2:$bpacc 3:$puacc\n";
	return ($acc, $bpacc, $puacc);
}

##################################################################################
# Calculates the base-pair accuracy of the given structure within a reference 
# sequence, although the structure can also be across the entire length of the 
# reference sequence. The accuracy is calculated according to the given dotplot, 
# i.e. base-pair probabilities or averaged probabilities for window-based approaches. The unpaired probabilities 
# The base-pair accuracy is actually the mean base-pair score for the given
# structure, i.e. sum of base-pair (pseudo-)probabilities divided by the number
# of base-pairs in the structure. The notation of the structures must be given
# in either alphabet defined in %STR_FORMATS. The default is the dotbracket format
# in $DOTBRACKTET.
# Also there is the option of calculating ranked accuracies to make different
# methods comparable. Here all base-pair probabilities within the dotplot are
# ranked after being rounded to x decimal places. Any identical values receive 
# identical ranks. The rounding to x.d.p. ensures that outliers do not heavily
# influence the average rank score. I.e. if some base-pairs are not predicted.
# but most are, then the rank can not be extremely high, but only up to 1o^x.
# NOTE: ALL INDICES COUNT FROM 1 TO N, except of course within arrays. This is due
#		to the counting in the Vienna package files.
# NOTE: This is only done for the base-pairs of the structure. Unpaired bases are
# ignored.
#
# Input: 
#		struct		The structure string in the one of the %STR_FORMATS alphabets.
#		str_f		(Optional) Format of structure in struct; one of %STR_FORMATS. 
#					Default is $DOTBRACKET, the dot-bracket format.
#		start		This is the starting position of the structure in the sequence.
#					(counting 1 to n)
#		bp_prob		A hash of ALL base-pair probabilities for the reference sequence,
#					where the key is "i:j" and i < j and they correspond to the given 
#					sequence.
#		seq			(Optional) The original/reference sequence. 
#					The structure is checked to see if it feasible given
#					this sequence (i.e. can the base-pairs form and are the
#					indices within the correct range?). 
#		seq_len		(Compulsory) Give the reference sequence length 
#					Also the indices of
#					the PU values are checked to be within the sequence boundaries.
#		rank		(Optional) 
#					A value x for the rounding of base-pair probabilities to x.d.p.
#					If this option is given it means ranked accuracies are calculated,
#					i.e. average base-pair ranks instead of using the raw value.
#		constrainToL This is the L parameter for which the sequence structure was
#					predicted. If this parameter is greater 0, then only base-pairs
#					are averaged that have a max. base-pair span of L, so 
#					calculating the accuracy of the "possible" structure. 
# Output:
#		The accuracy score or the ranked accuracy depending on whether rank is given
# 		and then the same for the possible structure given by L.
#		NOTE: if constrainToL == 0, the same value is given twice.
#
##################################################################################
sub BPaccuracy{
	my ($struct, $str_f, $start, $bp_prob_href, $seq, $seq_len, $rank, 
		$constrainToL) = @_;
	my $FUNCTION = "BPaccuracy in Structure.pm";
	
	# check input
	die ("INPUT ERROR in $FUNCTION:\nWhere does the structure start in you sequence?\n") unless ($start);
	die ("INPUT ERROR in $FUNCTION:\nWhere is the input structure?\n") unless ($struct);
	die ("INPUT ERROR in $FUNCTION:\nWhere are the base-pair probabilites?\n") unless ($bp_prob_href);
	if($seq){
		$seq_len = length($seq) unless ($seq_len);
	} else {
		die ("INPUT ERROR in $FUNCTION:\nThe reference sequence length ".
			"was not given!\n") unless ($seq_len);
	}
	
	## check that the structure can still fit into the sequence, wrt its start pos
	my $str_len = length($struct);
	my $str_end = $start + $str_len - 1;
	if($str_end > $seq_len){
		die ("INPUT ERROR in $FUNCTION: the structure cannot fit into the given".
			"sequence w.r.t the given starting position!\n");
	}
	
	# accuracy score - init
	my $bpacc = 0;
	my $bpaccL = 0;
	
	# get base-pairs from structure	counting 1 to n
	$str_f = $CONFIG{STR_FORMATS}->{DOTBRACKET} unless ($str_f);
	my $bp_href = structure2basepairs($struct, $start, $str_f);
	
	# check structure
	check_bps_wrt_seq($seq, $bp_href) if($seq);
	
	## cannot have no base-pairs, this wouldn't make sense
	my @bp_str = keys(%{$bp_href});
	die "ERROR in $FUNCTION: this is not a valid structure! $struct\n" 
		if(@bp_str < 2);
		
	## if rank is given, round the base-pair probabilities and calculate
	## ranks for each base-pair.
	my $highest_rank = 10**$rank;
	if($rank){
		## convert probabilities (rounded to x.d.p.) in $bp_prob_href to ranks
		$highest_rank = rank_bp_href($bp_prob_href, $rank);
	} 
	
	my $num_bp_strL = 0;
	my $num_bp = 0;
	## sum of base-pair probabilities/ranks
	foreach my $bp (@bp_str){
		++$num_bp;
		if(defined($bp_prob_href->{$bp})){
			$bpacc += $bp_prob_href->{$bp};
			if($constrainToL){
				my @bppos = split(":", $bp);
				if($bppos[1]-$bppos[0]+1 <= $constrainToL){
					++ $num_bp_strL;
					$bpaccL += $bp_prob_href->{$bp};
				}
			}
		## this is the case when the probability is 0!
		} else {
			if($rank){
				$bpacc += $highest_rank;
				
					if($constrainToL){
						my @bppos = split(":", $bp);
						if($bppos[1]-$bppos[0]+1 <= $constrainToL){
							++ $num_bp_strL;
							$bpaccL += $highest_rank;
						}
					}
			} else {
				if($constrainToL){
					my @bppos = split(":", $bp);
					if($bppos[1]-$bppos[0]+1 <= $constrainToL){
						++ $num_bp_strL;
					}
				}
			}
		}
	}
	# normalised (averaged) accuracy or ranked accuracy and the constrained accuracy
	my @res;
	if($constrainToL){
		if($num_bp_strL == 0){
			@res = ($bpacc/$num_bp, "NA");
		} else{
			@res = ($bpacc/$num_bp, $bpaccL/$num_bp_strL);
		}
	} else {
		@res = ($bpacc/$num_bp, $bpacc/$num_bp);		
	}
	
	## return 
	return @res;
}

##################################################################################
# Ranks the base-pair probabilities in the given hash reference. First sorts
# the probabilities so that the highest probability gets the highest rank.
# To avoid having very many ranks you can round the probabilities to x.d.p., 
# which means there are at most 10^x ranks. Also equal values get the
# same rank. Replaces the base-pair hash with the actual ranks instead of the
# probability.
#
# INPUT:
#		bp_prob_href	(Compulsory)
#						The base-pair probabilities for a sequence 
#						in a hash reference. Keys in format "i:j" and i<j.
#		round_xdp		(Optional)		
#						If given, the probabilities are rounded to the number
#						of decimal places of this parameter.
# OUTPUT:
#		None. The ranks are written into the input parameter base-pair hash
#		reference and overwrite the original base-pair probabilities.
##################################################################################
sub rank_bp_href{
	my ($bp_prob_href, $round_xdp) = @_;
	my $FUNCTION = "Structure::rank_bp_href";
	
	## check input variables
	die "INPUT ERROR in $FUNCTION:\nThere is no base-pair probability hash reference!"
		unless ($bp_prob_href);
	
	## hash to store unique (ranked) probabilities
	my %unique_probs = ();
	my @bps = keys (%{$bp_prob_href});
	
	## add all probabilities to the unique hash (set)
	foreach my $bp (@bps){
		if($round_xdp){
			my $rounded = StructureLibrary::Tools::round_Xdp($bp_prob_href->{$bp},
				$round_xdp);
				## change value in hash to rounded value
				$bp_prob_href->{$bp} = $rounded;
			$unique_probs{$rounded} = 0;
		} else {
			$unique_probs{$bp_prob_href->{$bp}} = 0;
		}
		
	}	
	
	## sort keys of hash, i.e. sort set
	my @sorted_probs = sort {$b<=>$a} (keys(%unique_probs));
	
	## map sorted array to ranks and have these as the value of the %unique_probs 
	for (my $i = 0 ; $i < @sorted_probs ; $i++){
		$unique_probs{$sorted_probs[$i]} = $i + 1;
	}
	
	## fill ranks into base-pair hash
	foreach my $bp (@bps){
		$bp_prob_href->{$bp} = $unique_probs{$bp_prob_href->{$bp}};
	}
	## this is the highest rank for the base-pair probability of 0.
	return @sorted_probs +1
}



##################################################################################
# Checks the base-pairs in the keys of a hash to see whether they are allowed
# according to the given sequence. Allowed base-pairs are A:U, U:A, G:C, C:G,
# G:U, and U:G. If a forbidden base-pair pairs, then an error message is reported
# and the program quits, otherwise 1 is returned.
# NOTE: also works on DNA sequences.
#
# INPUT:
#		seq			The template sequence.
#		bp_href		A hash reference containing base-pairs as keys in the format "i:j"
# OUTPUT:
#		1 if the base-pairs of the structure can form, given the sequence.
##################################################################################
sub check_bps_wrt_seq{
	my ($seq, $bp_href) = @_;
	my $FUNCTION = "check_bps_wrt_seq in Structure.pm";
	
	($seq) or die ("INPUT ERROR in $FUNCTION:\nThe sequence input seq is empty!\n");
	($bp_href) or die ("INPUT ERROR in $FUNCTION:\nThe structure hash input bp_href is empty!\n");
	
	my @seq_a = split("", uc($seq));
	
	# check each base-pair and die with an error message if it is not allowed
	# according to the given sequence
	my ($i, $j) = ();
	my ($base_i, $base_j) = "";
	foreach (keys (%{$bp_href})){
		if (($i, $j) = split(":", $_)){
			
			die ("ERROR in $FUNCTION:\n $i is out of bounds!\n") unless($base_i = $seq_a[$i-1]);
			die ("ERROR in $FUNCTION:\n $i is out of bounds!\n") unless($base_j = $seq_a[$j-1]);
			
			
			die("ERROR in $FUNCTION:\nThis base-pair is not allowed: $_!\n") 
					if($base_i eq "A" && !($base_j eq "U" || $base_j eq "T"));
			die("ERROR in $FUNCTION:\nThis base-pair is not allowed: $_!\n") 
					if(($base_i eq "U" || $base_i eq "T") && !($base_j eq "A" || $base_j eq "G"));
			die("ERROR in $FUNCTION:\nThis base-pair is not allowed: $_!\n") 
					if($base_i eq "G" && !($base_j eq "C" || $base_j eq "U"));
			die("ERROR in $FUNCTION:\nThis base-pair is not allowed: $_!\n") 
				if($base_i eq "C" && !($base_j eq "G"));
		} else {
			die("ERROR in $FUNCTION:\nThe base-pair keys in hash have the wrong format $_\n");
		}
	}
	#passed all checks
	return 1;
}

##################################################################################
# Removes base-pairs that are not compatible with the sequence.
#
# INPUT:
#		seq			The template sequence.
#		bp_href		A hash reference containing base-pairs as keys in the format "i:j"
# OUTPUT:
#		The base pairs in a hash-reference (same format as bp_href)
##################################################################################
sub remove_non_compatible_bps_wrt_seq{
	my ($seq, $bp_href) = @_;
	my $FUNCTION = "remove_non_compatible_bps_wrt_seq in Structure.pm";
	
	($seq) or die ("INPUT ERROR in $FUNCTION:\nThe sequence input seq is empty!\n");
	($bp_href) or die ("INPUT ERROR in $FUNCTION:\nThe structure hash input bp_href is empty!\n");
	
	my @seq_a = split("", uc($seq));	
	
	# check each base-pair and die with an error message if it is not allowed
	# according to the given sequence
	my ($i, $j) = ();
	my ($base_i, $base_j) = "";
	foreach (keys (%{$bp_href})){
		if (($i, $j) = split(":", $_)){
			die ("ERROR in $FUNCTION:\n $i is out of bounds!\n") unless($base_i = $seq_a[$i-1]);
			die ("ERROR in $FUNCTION:\n $i is out of bounds!\n") unless($base_j = $seq_a[$j-1]);
			
			delete($bp_href->{$_}) 
				if($base_i eq "A" && !($base_j eq "U" || $base_j eq "T"));
			delete($bp_href->{$_}) 
					if(($base_i eq "U" || $base_i eq "T") && !($base_j eq "A" || $base_j eq "G"));
			delete($bp_href->{$_}) 
					if($base_i eq "G" && !($base_j eq "C" || $base_j eq "U"));
			delete($bp_href->{$_})
				if($base_i eq "C" && !($base_j eq "G"));
			delete($bp_href->{$_})
				if($base_i eq "." || $base_j eq ".");
		} else {
			die("ERROR in $FUNCTION:\nThe base-pair keys in hash have the wrong format $_\n");
		}
	}
	#passed all checks
	return 1;
}

##################################################################################
# Parses the given PU file and returns all unpaired probabilities for ONE u, which
# is the stretch of bases that is unpaired. The file must be in the following
# format: an arbitrary number of header lines, that is defined by $num_header,
# then a line of probabilities for each index from 1 to n. The line has first the
# index number then the value for u=1, u=2, u=3,...u=max_u, separated by tabs.
# The index position represents the right border of an interval. E.g. for u=5 and
# and index of 10, the probability in the 6th column represents the probability that
# bases 6-10 are unpaired (a region of 5 bases).
# NOTE: with start and stop you can specify a region in which to collect the PU
# values.
#
# INPUT: 
#		PUfile		The file containing the unpaired probabilites.
#		u			The size of the unpaired region
#		num_header	The number of header lines in the PUfile
#		start		(Optional) Give a starting index (1-n) from which to fetch PUs
#		stop		(Optional) Give an ending index (1-n) until which to fetch PUs
# OUTPUT:
#		An array ref for the probability of a stretch of u bases to be unpaired for
#		entire sequence or specified region, counting 0 to n-1. 
#		Values can equal "NA" when no probability exists,
#		because no large enough region exists at the beginning of the sequence.
#		For example the single-based unpaired probability for base i, you will find
#		at position (i-1) of the output array. 
##################################################################################
sub parse_PUfile_get1u{
	my ($PUfile, $u, $num_header, $start, $stop) = @_;
	my $FUNCTION = "parse_PUfile_get1u in Structure.pm";
		
	open(IN_HANDLE, "<$PUfile") || die "ERROR in $FUNCTION:\n Couldn't open file $PUfile/n";
	
	$start = 1 unless ($start);
	
	my @pu = ();
	my $header_count = 0;
	my @row = ();
	my $index = 0;
	while(my $line = <IN_HANDLE>){
		chomp($line);
		++$header_count;
		
		# ignore header lines
		if($header_count > $num_header){
			@row = split("\t", $line);
			$index = $row[0];
			if(!$stop && $index >= $start){
				push(@pu, $row[$u]);
			} elsif($index >= $start && $index <= $stop){
				push(@pu, $row[$u]);
			}
		}
	}
	close(IN_HANDLE);
	die "WARNING in $FUNCTION: There are no accessibility values in the file".
		" $PUfile!\n" unless (@pu);
	return @pu;	
}

##################################################################################
# Parses the given PU file and returns all unpaired probabilities for GIVEN u's, 
# where u is the stretch of bases that is unpaired. The file must be in the following
# format: an arbitrary number of header lines, that is defined by $num_header,
# then a line of probabilities for each index from 1 to n. The line has first the
# index number then the value for u=1, u=2, u=3,...u=max_u, separated by tabs.
# The index position represents the right border of an interval. E.g. for u=5 and
# and index of 10, the probability in the 6th column represents the probability that
# bases 6-10 are unpaired (a region of 5 bases).
# NOTE: with start and stop you can specify a region in which to collect the PU
# values.
#
# INPUT: 
#		PUfile		The file containing the unpaired probabilites.
#		us_afref	A list of u's, sizes of unpaired regions to fetch
#		num_header	The number of header lines in the PUfile
#		start		(Optional) Give a starting index (1-n) from which to fetch PUs
#		stop		(Optional) Give an ending index (1-n) until which to fetch PUs
# OUTPUT:
#		A hash reference where they key is equal one of the given u's (size of 
#		unpaired region) and the value is an array reference for all unpaired
#		probabilities for that u across the entire sequence, counting from 
#		0 to n-1, where the index signifies the right border of the region.
#		For this reason, there are possible "NA" values at the beginning of the
#		array. So to get at base 50 for u=1 you need hashref-{1}->[49].
##################################################################################
sub parse_PUfile_getsomeu{
	my ($PUfile, $us_aref, $num_header, $start, $stop) = @_;
		my $FUNCTION = "parse_PUfile_get1u in Structure.pm";
		
	open(IN_HANDLE, "<$PUfile") || die "ERROR in $FUNCTION:\n Couldn't open file $PUfile/n";
	
	$start = 1 unless ($start);
	
	# initialise output hash
	my %pu = ();
	foreach (@{$us_aref}){
		my @u_row = ();
		$pu{$_} = \@u_row;
	}
	
	my @pu = ();
	my $header_count = 0;
	my @row = ();
	my $index = 0;
	while(my $line = <IN_HANDLE>){
		chomp($line);
		++$header_count;
		
		# ignore header lines
		if($header_count > $num_header){
			@row = split("\t", $line);
			$index = $row[0];
			if(!$stop && $index >= $start){
				foreach (@{$us_aref}){
					push(@{$pu{$_}}, $row[$_]);
				}
			} elsif($index >= $start && $index <= $stop){
				foreach (@{$us_aref}){
					push(@{$pu{$_}}, $row[$_]);
				}
			}
		}
	}
	close(IN_HANDLE);
	return \%pu;
}

##################################################################################
# Parses the given PU file and returns all unpaired probabilities for ALL u's, 
# where u is the stretch of bases that is unpaired. The file must be in the following
# format: an arbitrary number of header lines, that is defined by $num_header,
# then a line of probabilities for each index from 1 to n. The line has first the
# index number then the value for u=1, u=2, u=3,...u=max_u, separated by tabs.
# The index position represents the right border of an interval. E.g. for u=5 and
# and index of 10, the probability in the 6th column represents the probability that
# bases 6-10 are unpaired (a region of 5 bases).
# NOTE: with start and stop you can specify a region in which to collect the PU
# values.
#
# INPUT: 
#		PUfile		The file containing the unpaired probabilites.
#		u			The size of the unpaired region
#		num_header	The number of header lines in the PUfile
#		start		(Optional) Give a starting index (1-n) from which to fetch PUs
#		stop		(Optional) Give an ending index (1-n) until which to fetch PUs
# OUTPUT:
#		An array of arrays as a reference of unpaired probabilities for all regions,
#		i.e. u's. The first array indicates the position in the sequence, 
#		and the second array at index i contains probabilites for all
#		u's in the file, counting form 0 to n-1. So to get at base 50 for u=1 
#		you need arrayref->[49]->[0].
##################################################################################
sub parse_PUfile_getallu{
	my ($PUfile, $num_header, $start, $stop) = @_;
		my $FUNCTION = "parse_PUfile_get1u in Structure.pm";
		
	open(IN_HANDLE, "<$PUfile") || die "ERROR in $FUNCTION:\n Couldn't open file $PUfile/n";
	
	$start = 1 unless ($start);
	
	my @pu = ();
	my $header_count = 0;
	my @row = ();
	my $index = 0;
	while(my $line = <IN_HANDLE>){
		chomp($line);
		++$header_count;
		
		# ignore header lines
		if($header_count > $num_header){
			@row = split("\t", $line);
			$index = $row[0];
			if(!$stop && $index >= $start){
				shift(@row);
				push(@pu, \@row);
			} elsif($index >= $start && $index <= $stop){
				shift(@row);
				push(@pu, \@row);
			}
		}
	}
	close(IN_HANDLE);
	return \@pu;
}

##################################################################################
# INITIAL CHECK PASSED, COUNTS CORRECTLY AND 1 to n
# This method parses a structure in one of the formats given in %STR_FORMATS 
# notation and converts it to a list of base-pairs and returns this. 
# NOTE: you must given the symbol alphabet
# Only round brackets are considered as being base-pairs, the other symbols will 
# be checked to see if they fit to the format, but otherwise ignored.
# NOTE: The base-pair indices are given from 1 to n, not 0 to n-1 as this is
# how the output of the Vienna package is given as well.
# INPUT: 
#		db			Structure in dot-bracket notation.
#		start		(Optional) The starting index of the dot-bracket notation 
#					in the original sequence. Given as 1 to n counting.
#		str-f 		The structure format, choose one defined in global
#					variable %STR_FORMATS.
#		check		(Optional) First check whether the structure fits to the given format.
# OUTPUT: 
# 		A hash reference with the base-pairs as the keys in the form of "i:j",
#		i<j and values equals 1, so you can check whether a base-pair exists in the 
#		structure.
################################################################################
sub structure2basepairs{
	my ($struct, $start, $str_f, $check) = @_;
	my $FUNCTION = "structure2basepairs in Structure.pm";
	
	# structure format
	my @symba = ();
	check_str_format($struct,$str_f,$FUNCTION) if ($check);
	@symba = split("", $str_f);
	
	#start
	$start = 1 unless($start);
	
	# init base-pair hash
	my %bp = ();
	
	# parse structure
	my @struct_a = split("", $struct);
	my $n = @struct_a;
	my (@open, @close) = ();
		
	# parse dot-bracket notation
	for (my $i = 0; $i < $n; $i++){
		if($struct_a[$i] eq "("){
			push(@open, ($i+$start));
		} elsif($struct_a[$i] eq ")"){
			if(@open){
				$bp{(pop(@open)).":".($i+$start)} = 1;
			} else {
				die("ERROR in $FUNCTION: The given structure has the wrong format, i.e. \n".
				"there are an uneven number of opening and closing brackets: $struct!\n");
			}
		} 
	}
	die("ERROR in $FUNCTION: The given structure has the wrong format, i.e. \n".
			"there are an uneven number of opening and closing brackets: $struct!\n") if(@open);

	# return base-pairs
	return \%bp;
}


##################################################################################
# INITIAL CHECK PASSED, COUNTS CORRECTLY AND 1 to n
# This method gets the base-pairs in a hash reference and changes this to 
# the dot-bracket format.
#
# INPUT: 
#		bp_href		The hash reference of base-pairs in the form of 
#		len			The length of the structure (this is important to know
#					how many unpaired positions are at the end!)
# OUTPUT: 
# 		A hash reference with the base-pairs as the keys in the form of "i:j",
#		i<j and empty values, so you can check whether a base-pair exists in the 
#		structure.
################################################################################
sub basepairs2dotbracket{
	my ($bp_href, $len) = @_;
	my $FUNCTION = "structure2basepairs in Structure.pm";
	
	die "ERROR in $FUNCTION: Input parameters are missing!\n" unless ($bp_href, $len);
	
	## initiate stucture array for dots and brackets, initiate all 
	## elements with a dot
	my @struct_db = ();
	foreach (0..($len-1)){
		push(@struct_db, ".");
	}
	
	## add base pairs, i.e. brackets
	my ($i, $j) = 0;
	foreach (keys (%{$bp_href})){
		if (($i, $j) = split(":", $_)){
			if($i < $j){
				$struct_db[($i - 1)] = "(";
				$struct_db[($j - 1)] = ")";
			} elsif ($i > $j) {
				print STDERR "WARNING in $FUNCTION: The base-pair hash is in the".
					" wrong format as i is greater than j in the key. Swapping i ".
					"and j.\n";
				$struct_db[($j - 1)] = "(";
				$struct_db[($i - 1)] = ")";
			}
		} else {
			die("ERROR in $FUNCTION:\nThe base-pair keys in hash have the wrong format $_\n");
		}
	}
	
	return (join('',@struct_db));
}


##################################################################################
# This method just checks whether the input structure fits to the given format
# and whether the format is one of the known ones defined in %STR_FORMATS.
# It could be that the structure format is a subset of one of the known formats,
# this is also checked and the correct (full) format is given.
# Input:
#		struct		(Optional) The structure string, if not given, it is not checked
#					if the structure string contains the correct symbols.
#		str_f		The structure format
#		function	The name of the function from where you are checking
# Output:
#		The full format, i.e. symbols if the check is passed.
##################################################################################
sub check_str_format{
	my ($struct, $str_f, $function) = @_;
	
	my $format = 0;
	# check that format is allowed
	my $yes = 0;
	foreach (keys %{$CONFIG{STR_FORMATS}}){
		$format = $CONFIG{STR_FORMATS}->{$_} if ($str_f eq $CONFIG{STR_FORMATS}->{$_});
	}
	
	# maybe it is in the wrong order, or only a subset of the format is used (also ok)
	unless ($format){
		my @str_f_a = split ("", $str_f);
		my $n = @str_f_a; 
		foreach (keys %{$CONFIG{STR_FORMATS}}){
			my $check = 1;
			for(my $i = 0 ; $i < $n ; $i++){
				if(index($CONFIG{STR_FORMATS}->{$_},$str_f_a[$i]) == -1 ){
					$check = 0;
				} 
			}
			$format = $CONFIG{STR_FORMATS}->{$_} if ($check);
		}
	}
	
	($format) or die "ERROR in $function: The given structure format, $str_f, is unknown!\n";
	
	# check that structure consists of the given symbols
	my $yes = 0;
	if($struct){
		$yes = 0;
		foreach (keys %{$CONFIG{STR_FORMATS}}){
			if ($struct =~ /^[$format]+$/){
			}
			$yes = 1 if ($struct =~ /^[$format]+$/);
		}
		($yes) or die "ERROR in $function: The given structure does not match the structure format!\n".
						"Format: $str_f\nStructure:$struct\n";
	}
					
	return $format;
}


##################################################################################
# This method is given two fasta files with the identical sequence (not checked).
# The second file also contains a constraint folding of the sequence. The sequence
# is then folded by RNAfold first normally and then with the constraint. Then the
# energies are parsed and the energy difference (ED) calculated and the probability
# of the constraint structure is also given. The results are in a hash reference
# as defined below.
# Input:
#		fasta_ori 	The fasta file just including the header and the sequence
#					where the sequence MUST be in one line! This is not checked!!
#		fasta_const	The fasta file with the sequence, but after the sequence
#					the RNAfold constraint sequence is given. NOTE: the constraint
#					sequence must be compatible with RNAfold and the same length
#					as the sequence. This is also not checked!!
#		rnafold		This is the path to the RNAfold installation. It is best to
#					use the complete path and not just RNAfold.
#		params		(Optional) This is a list of parameters that are compatible 
#					with RNAfold. If it is not given, the default parameters
#					are used.
#		dir			(Optional) A location can be given for writing the RNAfold output. 
#					This can be useful if there are writing issues due to the cluster
#					or there is a specific location to which the results should be saved.
#					If dir is not given, then the file paths are relative to the 
#					current directory.
#		nops		(Optional) No ps files are generated when this option if defined, 
#					if 0 or undefined, then just the ensemble energies are calculated.
#		silent		(Optional) Suppress all prints, even warnings if defined.
#		prevres_href (Optional) Maybe you have already called RNAfold for this 
#					sequence and you just need to calculate a new constraint energy.
#					Then CONSTRAINT_E, ED, and CONSTRAINT_P are recalculated and
#					updated in the hash.
# Output:
#		A hash reference with the output information from RNAfold as follows
#		Key			Value
#		MFE			The MFE structure without constraint
#		MFE_CONSTRAINT	The MFE structure with constraint
#		ENSEMBLE_E	The ensemble energy of the sequence as written to STDOUT by RNAfold
#		CONSTRAINT_E The ensemble energy of all structures that fit to the constraint
#		MFE_E		The energy of the MFE structure 
#		CONSTRAINT_MFE_E The energy of the constraint MFE structure
#		ED			The energy difference between all structures and constraint structures
#		CONSTRAINT_P The probability of the constraint structures
#		CONSTRAINT_MFE_P The probability of the constraint MFE structure
#		DP			The ps file destination of the dotplot without constraint
#		SS			The ps file destination of the structure plot without constraint
##################################################################################
sub ensemble_energies_and_probability_with_constraint{
	my ($fasta_ori, $fasta_const, $params, $dir, $nops, $silent, $prevres_href) = @_;
	my $FUNCTION = "ensemble_energies in Structure.pm";
	
	my $RT = $CONFIG{RT};
	
	# call RNAfold and get results 
	my $res_all_href = $prevres_href; # hash reference for all output results
	$res_all_href->{CONSTRAINT_E} = "NA";
	$res_all_href->{ED} = "NA";
	$res_all_href->{CONSTRAINT_P} = "NA";
	($prevres_href) or $res_all_href = call_RNAfold_parse_output($fasta_ori, $params, $dir, $nops, $silent);
	$params .= " -C";
	my $res_constraint_href = call_RNAfold_parse_output( $fasta_const, $params, $dir, $nops, $silent);
	
	# add structure
	$res_all_href->{CONSTRAINT_MFE} = $res_constraint_href->{MFE};
	
	# add constraint ensemble energy to the results hash
	$res_all_href->{CONSTRAINT_E} = $res_constraint_href->{ENSEMBLE_E};
	
	# add constraint mfe energy to the results hash
	$res_all_href->{CONSTRAINT_MFE_E} = $res_constraint_href->{MFE_E};
	
	# calculate ED
	my $ed = $res_all_href->{ENSEMBLE_E} - $res_all_href->{CONSTRAINT_E};
	$res_all_href->{ED} = $ed;
	
	# calculate constraint MFE ED
	my $mfe_ed = $res_all_href->{ENSEMBLE_E} - $res_all_href->{CONSTRAINT_MFE_E};
	$res_all_href->{CONSTRAINT_MFE_P} = exp($mfe_ed/$RT);
	
	# calculate probability
	$res_all_href->{CONSTRAINT_P} = exp($ed/$RT);
	
	return $res_all_href;
}

##################################################################################
# Calls the given installation of RNAfold with the given options or the 
# default options which are: RNAfold -p -d2 -noLP
# Returned is a hash reference with the relevant information as defined below.
# Input:
#		fasta 		The fasta file including the sequence in ONE LINE and 
#					if constraint folding is used, then this constraint sequence
#					is on the line following the sequence.
#		params		(Optional) This is a list of parameters that are compatible 
#					with RNAfold. If it is not given, the default parameters
#					are used.
#		dir			(Optional) A location can be given for writing the RNAfold output. 
#					This can be useful if there are writing issues due to the cluster
#					or there is a specific location to which the results should be saved.
#					If dir is not given, then the file paths are relative to the 
#					current directory.
#		nops		(Optional) No ps files are generated when this option if defined, 
#					if 0 or undefined, then just the ensemble energies are calculated.
#		noprint		(Optional) Suppress all prints, even warnings if defined.
# Output:
#		A hash reference with the output information from RNAfold as follows
#		KEY			VALUE
#		ENSEMBLE_E	The ensemble energy as written to STDOUT by RNAfold
#		MFE_E		The energy of the MFE structure 
#		MFE			The MFE structure
#		DP			The ps file destination of the dotplot
#		SS			The ps file destination of the structure plot
##################################################################################
sub call_RNAfold_parse_output{
	my ($fasta, $params, $dir, $nops, $silent) = @_;	
	my $FUNCTION = "call_RNAfold_parse_output in Structure.pm";
	
	 my $fasta = abs_path($fasta);
	
	# results hash
	my %results = ();
	
	# parameters are not checked, this is done implicitly further down 
	
	# tmp directory
	my $tmpdir = "";
	if($dir){
		$dir .= "/";
		$tmpdir = $dir."tmp_rnafold/"
	} else {
		$dir = ".";
		$tmpdir = $dir."/tmp_rnafold/";
		
	}
	(-e $tmpdir) or system ("mkdir $tmpdir");
	
	
	# RNAfold call according to the input parameters
	my $RNAfoldcall = "cd $tmpdir ; $CONFIG{RNAFOLD} ";
	if($params){
		$RNAfoldcall .= $params;
	} else {
		$params = "-p -d2 -noLP";
		print STDERR "WARNING: using default parameters for RNAfold $params\n" unless($silent);
		$RNAfoldcall .= $params;
	} 
	if ($nops){
		$RNAfoldcall .= " -noPS "; 
	}
	$RNAfoldcall .= " < $fasta";
		
	# call RNAfold	
	my @RNAfold_results = readpipe($RNAfoldcall);
	
	# only deal with output files if they are wanted
	if(!defined($nops)){
		# declare output file names
		my $dp_file = "";
		my $ss_file = "";
	
		# read tmp directory
		opendir(IMD, "$tmpdir") || die("ERROR: Cannot open directory $tmpdir");
			my @files= readdir(IMD);
			foreach my $file (@files){
				if($file =~ /(.+_dp\.ps)$/){
					$dp_file = $tmpdir.$1;
				} elsif ($file =~ /(.+_ss\.ps)$/){
				$ss_file = $tmpdir.$1;
				} else {
					print "WARNING: what is this file???? $file\n" unless(
						$file eq "." || $file eq ".." );
				}
			}
		closedir(IMD);
		
		# move result files and delete tmpdir
		system("mv $dp_file $dir");
		die ("ERROR in $FUNCTION: could not move $dp_file!\n") 
			unless($dir.$dp_file);
		$dp_file = $dir.$dp_file;
		system("mv $dp_file $dir");	
		die ("ERROR in $FUNCTION: could not move $ss_file!\n") 
			unless(-e $dir.$ss_file);
		$ss_file = $dir.$ss_file;
	
		# delete temporary directory
		system("rm -R $tmpdir");
		die ("ERROR in $FUNCTION: could not delete $tmpdir!\n") 
			if(-e $tmpdir);


		# save files in hash
		$results{DP} = $dp_file;
		$results{SS} = $ss_file;
	}
	
	# get energies
	my ($e_all, $e_mfe, $mfe) = "";
	foreach (@RNAfold_results) {
		chomp($_);
	    # get line with ensemble-energy 
	    if ($_ =~ /.+\[\s*(\-\d+\.\d+)\]$/) {
	       	 	$e_all = $1;
	    # get mfe and mfe energy
	    } elsif ($_ =~ /([\(\)\.]+)\s+\(\s*(\-?\d+\.\d+)\)$/){
				$mfe = $1;
	    		$e_mfe = $2;
	    # the mfe is completely unstructured
	    } elsif ($_ =~ /(\.+)\s+\(\s*(\-?\d+\.\d+)\)$/){
	    		$mfe = $1;
	    		$e_mfe = $2;
	    }
	}
	
	# If the call above did not include the energies, then something went wrong!
	if(!$e_all || !$e_mfe || !$mfe){
		my $RNAout = "";
		foreach (@RNAfold_results){
			$RNAout .= $_."\n";
		}
		print "ENSEMBLE ENERGY: $e_all\n";
		print "MFE ENERGY: $e_mfe\n";
		print "MFE: $mfe\n";
		die "ERROR in $FUNCTION: see RNAfold output below!\n $RNAout\n";
	}
	
	# add energies to hash and return it
	$results{ENSEMBLE_E} = $e_all;
	$results{MFE_E} = $e_mfe;
	$results{MFE} = $mfe;
	return \%results;
}


##################################################################################
# Calls the installation of RNAfold with the given options or the 
# default options (see RNAtoolsConfig for options). The results can be written
# to a given destination directory also with a given prefix before the normal
# endings _dp.ps for the dotplot and _ss.ps for the MFE structure file. 
# Returns an array of the resulting ps file names, with complete path.
# NOTE: there should only be one sequence in the fasta file and in one line!
# Input:
#		fasta 		The fasta file including the sequence in ONE LINE and 
#					if constraint folding is used, then this constraint sequence
#					is on the line following the sequence.
#		params		(Optional) This is a list of parameters that are compatible 
#					with RNAfold. If it is not given, the default parameters
#					are used. 
#		dir			(Optional) A location can be given for writing the RNAfold 
#					output. This can be useful if there are writing issues due 
#					to the cluster or there is a specific location to which the 
#					results should be saved. If dir is not given, then the file 
#					paths are relative to the current directory.
#		prefix		Give the output ps files this prefix before _dp.ps and _ss_ps.
# Output
#		An array of the resulting ps files with complete file path.
#		1. Doplot "_dp.ps"
#		2. MFE structure plot "_ss.ps"
##################################################################################
sub call_RNAfold{
	my ($fasta, $params, $dir, $prefix) = @_;	
	my $FUNCTION = "call_RNAfold in Structure.pm";
	
	return "" unless (-e $fasta);
	
	my $fasta = abs_path($fasta);
	
	# results array
	my @results = ();
	
	# parameters are not checked, this is done implicitly further down 
	
	# tmp directory
	my $tmpdir = "";
	if($dir){
		$dir .= "/";
		$tmpdir = $dir."tmp_rnafold/"
	} else {
		$dir = ".";
		$tmpdir = $dir."/tmp_rnafold/";
		
	}
	(-e $tmpdir) or system ("mkdir $tmpdir");
	
	
	# RNAfold call according to the input parameters
	my $RNAfoldcall = "cd $tmpdir ; $CONFIG{RNAFOLD} ";
	if($params){
		$RNAfoldcall .= $params;
	} else {
		$params = $CONFIG{RNAFOLDPARAMS};
		print STDERR "WARNING: using default parameters for RNAfold $params\n";
		$RNAfoldcall .= $params;
	} 
	$RNAfoldcall .= " < $fasta";
		
	# call RNAfold	
	my @RNAfold_results = readpipe($RNAfoldcall);

	# declare output file names
	my $dp_file = "";
	my $ss_file = "";
	
	# read tmp directory
	opendir(IMD, "$tmpdir") || die("ERROR: Cannot open directory $tmpdir");
		my @files= readdir(IMD);
		foreach my $file (@files){
			if($file =~ /(.+_dp\.ps)$/){
				$dp_file = $tmpdir.$1;
			} elsif ($file =~ /(.+_ss\.ps)$/){
				$ss_file = $tmpdir.$1;
			} else {
				print "WARNING: what is this file???? $file\n" unless(
					$file eq "." || $file eq ".." );
			}
		}
	closedir(IMD); 
	
	## if both files do not exist, there must have been
	## an error with RNAfold, so return -1, but do not die.
	return -1 unless(-e $dp_file && -e $ss_file);
	
	## remove automatic RNAfold name and replace with prefix
	## move to parent directory
	if($prefix){ # 
		if($dp_file =~ /_dp\.ps$/){
			system("mv $dp_file $dir/${prefix}_dp.ps");
			die "ERROR in $FUNCTION: could not move $dp_file!\n" unless
				(-e "$dir/${prefix}_dp.ps");
			$dp_file = "$dir/${prefix}_dp.ps";
		} else {
			die "ERROR in $FUNCTION: $dp_file did not have the format \\S+_dp\\.ps\n";
		}
		if($ss_file =~ /_ss\.ps$/){
			system("mv $ss_file $dir/${prefix}_ss.ps");
			die ("ERROR in $FUNCTION: could not move $ss_file!\n") 
				unless(-e "$dir/${prefix}_ss.ps");
			$ss_file = "$dir/${prefix}_ss.ps";
		} else {
			die "ERROR in $FUNCTION: $ss_file did not have the format \\S+_dp\\.ps\n";
		}
	## move result files to parent directory
	}  else {
		system("mv $dp_file $dir");
		die ("ERROR in $FUNCTION: could not move $dp_file!\n") 
			unless($dir.$dp_file);
		$dp_file = $dir.$dp_file;
		system("mv $ss_file $dir");	
		die ("ERROR in $FUNCTION: could not move $ss_file!\n") 
			unless(-e $dir.$ss_file);
		$ss_file = $dir.$ss_file;
	}
	
	# delete temporary directory
	system("rm -R $tmpdir");
	die ("ERROR in $FUNCTION: could not delete $tmpdir!\n") 
		if(-e $tmpdir);
	
	## return with full path
	return ($dp_file, $ss_file);
}

##################################################################################
# Calls the installation of RNAplfold with the given options or the 
# default options (see RNAtoolsConfig for options). The results can be written
# to a given destination directory also with a given prefix before the normal
# endings _dp.ps for the dotplot and _lunp for the accuracies file. 
# Returns an array of the resulting file name(s), with complete path.
# NOTE: there should only be one sequence in the fasta file and in one line!
# Input:
#		fasta 		The fasta file including the sequence in ONE LINE and 
#					if constraint folding is used, then this constraint sequence
#					is on the line following the sequence.
#		params		(Optional) This is a list of parameters that are compatible 
#					with RNAfold. If it is not given, the default parameters
#					are used. If the -u parameter is not given, then no
#					accessibilities are calculated, also not in the default.
#		dir			(Optional) A location can be given for writing the RNAfold 
#					output. This can be useful if there are writing issues due 
#					to the cluster or there is a specific location to which the 
#					results should be saved. If dir is not given, then the file 
#					paths are relative to the current directory.
#		prefix		Give the output ps files this prefix before _dp.ps and _lunp.
# Output
#		An array of the resulting dotplot and accuracy files with complete file path.
#		1. Doplot "_dp.ps"
#		2. Accuracies ".acc" or "_lunp" (default) and only if u > 0
##################################################################################
sub call_RNAplfold{
	my ($fasta, $params, $dir, $prefix) = @_;	
	my $FUNCTION = "call_RNAplfold in Structure.pm";
	
	return "" unless (-e $fasta);
	my $fasta = abs_path($fasta);
	
	# results array
	my @results = ();
	
	# parameters are not checked, this is done implicitly further down 
	
	# tmp directory
	my $tmpdir = "";
	if($dir){
		$dir .= "/";
		$tmpdir = $dir."tmp_rnaplfold/"
	} else {
		$dir = ".";
		$tmpdir = $dir."/tmp_rnaplfold/";
		
	}
	(-e $tmpdir) or system ("mkdir $tmpdir");
	
	
	# RNAplfold call according to the input parameters
	my $RNAplfoldcall = "cd $tmpdir ; $CONFIG{RNAPLFOLD} ";
	if($params){
		$RNAplfoldcall .= $params;
	} else {
		$params = $CONFIG{RNAPLFOLDPARAMS};
		print STDERR "WARNING in $FUNCTION:\n".
			"using default parameters for RNAplfold --> $params\n";
		$RNAplfoldcall .= $params;
	} 
	$RNAplfoldcall .= " < $fasta";
		
	# call RNAplfold	
	my @RNAplfold_results = readpipe($RNAplfoldcall);

	# declare output file names
	my $dp_file = "";
	my $acc_file = "";
	
	# read tmp directory
	opendir(IMD, "$tmpdir") || die("ERROR: Cannot open directory $tmpdir");
		my @files= readdir(IMD);
		foreach my $file (@files){
			if($file =~ /(.+_dp\.ps)$/){
				$dp_file = $tmpdir.$1;
			} elsif ($file =~ /(.+_lunp)$/){
				$acc_file = $tmpdir.$1;
			} else {
				print "WARNING: what is this file???? $file\n" unless(
					$file eq "." || $file eq ".." );
			}
		}
	closedir(IMD); 
	
	## if the dotplot file does not exist, there must have been
	## an error with RNAplfold, so return -1, but do not die.
	return -1 unless(-e $dp_file);
	## if accuracies should have been calculated, but no file exists, also return -1.
	if($params =~ /-u\s\d+/){
		return -1 unless (-e $acc_file);
	}
	
	## remove automatic RNAplfold name and replace with prefix
	## move to parent directory
	if($prefix){ # 
		if($dp_file =~ /_dp\.ps$/){
			system("mv $dp_file $dir/${prefix}.ps");
			die "ERROR in $FUNCTION:\ncould not move $dp_file!\n" unless
				(-e "$dir/${prefix}.ps");
			$dp_file = "$dir/${prefix}.ps";
		} else {
			die "ERROR in $FUNCTION:\n$dp_file did not have the format \\S+\\.ps\n";
		}
		if($acc_file && $acc_file =~ /_lunp$/){
			system("mv $acc_file $dir/${prefix}.acc");
			die ("ERROR in $FUNCTION: could not move $acc_file!\n") 
				unless(-e "$dir/${prefix}.acc");
			$acc_file = " $dir/${prefix}.acc";
		} else {
			die "ERROR in $FUNCTION: $acc_file did not have the format \\S+_lunp\n" 
				if ($acc_file);
		}
	## move result files to parent directory
	}  else {
		system("mv $dp_file $dir");
		die ("ERROR in $FUNCTION: could not move $dp_file!\n") 
			unless($dir.$dp_file);
		$dp_file = $dir.$dp_file;
		if($acc_file){
			system("mv $acc_file $dir");	
			die ("ERROR in $FUNCTION: could not move $acc_file!\n") 
				unless(-e $dir.$acc_file);
			$acc_file = $dir.$acc_file;
		}
	}
	
	# delete temporary directory
	system("rm -R $tmpdir");
	die ("ERROR in $FUNCTION: could not delete $tmpdir!\n") 
		if(-e $tmpdir);
	
	## return with full path
	return ($dp_file, $acc_file);
}

##################################################################################
# Calls the installation of the LocalFold (centred) program with the given 
# options or the default options (see RNAtoolsConfig for options). 
# The results can be written to a given destination directory also with
# a given prefix before the normal endings _dp.ps for the dotplot and _lunp for the accuracies file. 
# Returns an array of the resulting file name(s), with complete path.
# NOTE: there should only be one sequence in the fasta file and in one line!
# Input:
#		fasta 		The fasta file including the sequence in ONE LINE and 
#					if constraint folding is used, then this constraint sequence
#					is on the line following the sequence.
#		params		(Optional) This is a list of parameters that are compatible 
#					with RNAfold. If it is not given, the default parameters
#					are used. If the -u parameter is not given, then no
#					accessibilities are calculated, also not in the default.
#		dir			(Optional) A location can be given for writing the RNAfold 
#					output. This can be useful if there are writing issues due 
#					to the cluster or there is a specific location to which the 
#					results should be saved. If dir is not given, then the file 
#					paths are relative to the current directory.
#		prefix		Give the output ps files this prefix before _dp.ps and 
#					.acc or _lunp.
# Output
#		An array of the resulting dotplot and accuracy files with complete file path.
#		1. Doplot "_dp.ps"
#		2. Accuracies ".acc" or "_lunp" (default) and only if u > 0
##################################################################################
sub call_LocalFold{
	my ($fasta, $params, $dir, $prefix) = @_;	
	my $FUNCTION = "call_LocalFold in Structure.pm";
	
	return -1 unless (-e $fasta);
	my $fasta = abs_path($fasta);
	
	# results array
	my @results = ();
	
	# parameters are not checked, this is done implicitly further down 
	
	# tmp directory
	my $tmpdir = "";
	if($dir){
		$dir .= "/";
		$tmpdir = $dir."tmp_localfold/"
	} else {
		$dir = ".";
		$tmpdir = $dir."/tmp_localfold/";
		
	}
	(-e $tmpdir) or system ("mkdir $tmpdir");
	
	
	# RNAplfold call according to the input parameters
	my $localfoldcall = "cd $tmpdir ; $CONFIG{LOCALFOLD} ";
	if($params){
		$localfoldcall .= $params;
	} else {
		$params = $CONFIG{LOCALFOLDPARAMS};
		print STDERR "WARNING in $FUNCTION:\n".
			"using default parameters for LocalFold --> $params\n";
		$localfoldcall .= $params;
	} 
	$localfoldcall .= " -seqfile $fasta -noacc";
#	$localfoldcall .= " < $fasta"; ## TODO change when we have edited RNAplfold
		
	# call RNAplfold	
	system($localfoldcall);
#	my @localfold_results = readpipe($localfoldcall);

	# declare output file names
	my $dp_file = "";
	my $acc_file = "";
	
	# read tmp directory
	opendir(IMD, "$tmpdir") || die("ERROR: Cannot open directory $tmpdir");
		my @files= readdir(IMD);
		foreach my $file (@files){
			if($file =~ /(.+\.ps)$/){
				$dp_file = $tmpdir.$1;
			} elsif ($file =~ /(.+\.acc)$/){
#			} elsif ($file =~ /(.+_lunp)$/){ ## TODO change with edited RNAplfold
				$acc_file = $tmpdir.$1;
			} else {
				print "WARNING: what is this file???? $file\n" unless(
					$file eq "." || $file eq ".." );
			}
		}
	closedir(IMD); 
	
	## if the dotplot file does not exist, there must have been
	## an error with LocalFold, so return -1, but do not die.
	return -1 unless(-e $dp_file);
	## if accuracies should have been calculated, but no file exists, also return -1.
	if($params =~ /-u\s\d+/){
		return -1 unless (-e $acc_file);
	}
	
	## remove automatic RNAplfold name and replace with prefix
	## move to parent directory
	if($prefix){ # 
		if($dp_file =~ /\.ps$/){
			system("mv $dp_file $dir/${prefix}.ps");
			die "ERROR in $FUNCTION:\ncould not move $dp_file!\n" unless
				(-e "$dir/${prefix}.ps");
			$dp_file = "$dir/${prefix}.ps";
		} else {
			die "ERROR in $FUNCTION:\n$dp_file did not have the format \\S+\\.ps\n";
		}
		if($acc_file && $acc_file =~ /\.acc$/){
#		if($acc_file && $acc_file =~ /_lunp$/){ ## TODO change with edited RNAplfold
			system("mv $acc_file $dir/${prefix}.acc");
			die ("ERROR in $FUNCTION: could not move $acc_file!\n") 
				unless(-e "$dir/${prefix}.acc");
			$acc_file = " $dir/${prefix}.acc";
		} else {
			die "ERROR in $FUNCTION: $acc_file did not have the format \\S+\\.acc\n" 
				if ($acc_file);
#			die "ERROR in $FUNCTION: $acc_file did not have the format \\S+_lunp\n" 
#				if ($acc_file); ## TODO change with edited RNAplfold
		}
	## move result files to parent directory
	}  else {
		system("mv $dp_file $dir");
		die ("ERROR in $FUNCTION: could not move $dp_file!\n") 
			unless($dir.$dp_file);
		$dp_file = $dir.$dp_file;
		if($acc_file){
			system("mv $acc_file $dir");	
			die ("ERROR in $FUNCTION: could not move $acc_file!\n") 
				unless(-e $dir.$acc_file);
			$acc_file = $dir.$acc_file;
		}
	}
	
	# delete temporary directory
	system("rm -R $tmpdir");
	die ("ERROR in $FUNCTION: could not delete $tmpdir!\n") 
		if(-e $tmpdir);
	
	## return with full path
	return ($dp_file, $acc_file);
}


##################################################################################
# Reads a dotplot and returns the base-pair probabilities in the form of a
# hash reference, where the keys are the base pairs "i:j" (i < j) and the
# value is the probability.
# NOTE: in the dotplots, the indices count from 1 to n and not from 0 to n-1.
# Further note, that the values in the dotplots are the square root of the 
# probabilities, so here the values are squared.
#
# Input: 
#		dotplot	The dotplot file with base-pair probabilites as returned
#					by RNAfold or RNAplfold
# Output: 
#		A hash reference of the base-pair probabilites with "i:j" as the key,
#		i < j, and the probability as the value. 	
##################################################################################
sub parse_dotplot_return_bp_probs{
	my $dotplot = shift;
	my $FUNCTION = "parse_dotplot_return_bp_probs in Structure.pm";
	
	my %bp = ();
	
	open(IN_HANDLE, "<$dotplot") || die "ERROR in $FUNCTION:\n Couldn't open file $dotplot/n";
	
	while(my $line = <IN_HANDLE>){
		chomp($line);
		# base-pair probability line (works for colour converted dotplots too)
		# only for lines containting ubox
		if ($line =~ /ubox$/){
			if($line =~ /^(\d+)\s(\d+)\s(\S+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			# modified dotplots
			} elsif($line =~ /rgb.*\s(\d+)\s(\d+)\s(\S+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			}
		}
		
	}
	close(IN_HANDLE);
	return \%bp;
}

##################################################################################
# Reads a dotplot and returns the base-pair probabilities in the form of a
# hash reference, where the keys are the base pairs "i:j" (i < j) and the
# value is the probability.
# NOTE: in the dotplots, the indices count from 1 to n and not from 0 to n-1.
# Further note, that the values in the dotplots are the square root of the 
# probabilities, so here the values are squared.
#
# Input: 
#		dotplot	The dotplot file with base-pair probabilites as returned
#					by RNAfold or RNAplfold
# Output: 
#		An array ref including:
#	(1)	A hash reference of the base-pair probabilites with "i:j" as the key,
#		i < j, and the probability as the value. 
#	(2)	The sequence in the dotplot file as a string.
##################################################################################
sub parse_dotplot_return_bp_probs_and_sequence{
	my $dotplot = shift;
	my $FUNCTION = "parse_dotplot_return_bp_probs in Structure.pm";
	
	my %bp = ();
	my $seqline = ""; #lines in ps file including sequence
	my $inseq = 0;
	
	open(IN_HANDLE, "<$dotplot") || die "ERROR in $FUNCTION:\n Couldn't open file $dotplot/n";
	
	while(my $line = <IN_HANDLE>){
		chomp($line);
		
		# get sequence and convert it according to the given parameters
		if($line =~ /\/sequence/){
			$inseq = 1;
			$seqline = $line;
			} elsif($inseq){ # in sequence region
				if($line =~ /def/){ # end of sequence region
					$inseq = 0;
					$seqline .= $line;
				} else {
					$seqline .= $line;
			}
		
		# base-pair probability line (works for colour converted dotplots too)
		# only for lines containting ubox
		}elsif ($line =~ /ubox$/){
			if($line =~ /^(\d+)\s(\d+)\s(\S+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			# modified dotplots
			} elsif($line =~ /rgb.*\s(\d+)\s(\d+)\s(\S+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			}
		}
	}
	close(IN_HANDLE);
	
	# get rid of ps formatting
	$seqline =~ tr/\///d;
	$seqline =~ tr/\\//d;
	$seqline =~ tr/\(//d;
	$seqline =~ tr/\)//d;
	$seqline =~ s/sequence//;
	$seqline =~ s/def//;
	$seqline =~ tr/{//d;
	$seqline =~ tr/}//d;
	$seqline =~ tr/ //d;
	
	die "ERROR in $FUNCTION: The sequence is empty, either file is wrong,\n ".
		"or parsing failed in $dotplot\n" unless ($seqline);
	die "ERROR in $FUNCTION: There are no base-pair probabilities in $dotplot\n" unless (%bp);
	
	return (\%bp, $seqline);
}



##################################################################################
# Reads the base-pair probabilites as returned by "parse_dotplot_return_bp_probs"
# and adds the probabilites for each base to a single sum. This sum is the
# probability for the respective base to be PAIRED. A hash is returned,
# where the key is the index for the base and the value is the sum of all 
# probabilites of base-pairs in which that base is involved in, i.e.
# the probability for the base to be paired.
# NOTE: The numbering is from 1 to n and note that there are no entries for 
# positions that have 0 probability to be paired! You have to fill in these
# gaps yourself!
#
# Input:
#		bp_href		The hash reference to the base-pair probabilities
# Output:
#		A hash reference containing the base index as the key and the
# 		probability for that base to be paired as the value.
##################################################################################
sub prob_base_i_is_paired{
	my $bp_href = shift;
	my $FUNCTION = "prob_base_i_is_paired in Structure.pm";
	
	my %paired_prob = ();
	my ($i, $j) = 0;	
	foreach (keys (%{$bp_href})){
		if(my ($i,$j) = split(":", $_)){
			if(exists($paired_prob{$i})){
				$paired_prob{$i} += $bp_href->{$_};
			} else {
				$paired_prob{$i} = $bp_href->{$_};
			}
			if(exists($paired_prob{$j})){
				$paired_prob{$j} += $bp_href->{$_};
			} else {
				$paired_prob{$j} = $bp_href->{$_};
			}
		} else {
			die("ERROR in $FUNCTION:\nThe keys of the input hash are in the wrong format, e.g. $_\n");
		}
	}
#	foreach (keys (%paired_prob)){
#		print "paired prob of $_ is $paired_prob{$_}\n";
#	}
	
	return \%paired_prob;
}

##################################################################################
# Reads the fasta file including a header, a sequence and the corresponding 
# structure to the sequence and returns this as a hash reference.
# The file can include more than one item. The sequence can cover more than one
# line, but must contain [ATGCU]. The structure can also cover more than one
# line, but must contain the symbols given or "()." by default.
# NOTE: The method only checks to see if the sequence is the same length as the
# structure sequence, but not if the structure fits to the given sequence!
#
# INPUT:
#		file	The fasta file
#		str_f 	(Optional) The structure format. Default is $DOTBRACKET format.
# OUTPUT:
#		An array with 
#	(1)	A hash reference where the key is the item id and the value is another
#		hash reference with two items, SEQUENCE and STRUCTURE for the respective 
#		information. E.g. $href->{$id}->{SEQUENCE} or $href->{$id}->{STRUCTURE}.
#	(2)	An array reference including the ids in the order they are given in the
#		input file, $file. This information is necessary if you need the exact 
#		order, which is not given in the hash.
#
# NOTE: method function checked on 20.10.2010
##################################################################################
sub read_seq_structure_fasta_file{
	my($file, $str_f) = @_;
	my $FUNCTION = "read_seq_structure_fasta_file in Structure.pm";
	
	($str_f) or $str_f = $CONFIG{STR_FORMATS}->{DOTBRACKET};
	
	my $id 			= "";
	my $seqstring 	= "";
	my $structure 	= "";
	my %fasta 		= ();
	my $line 		= "";
	my @order 		= ();
	open(IN_HANDLE, "<$file") || die "ERROR in $FUNCTION:\nCouldn't open the following file in package Tool,".
									 " sub read_fasta_file: $file/n";
	
	while($line = <IN_HANDLE>){
		chomp($line);
		
		# header (can contain one space after > symbol)
		if($line =~ /^\>\s?(\S+)\s*/){
			if($id){
				my $lseq 	= 	length($seqstring);
				my $lstruc 	= 	length($structure);
				if($lseq   != 	$lstruc){
					die "ERROR IN $FUNCTION: The sequence length is not equal to the structure length in $file!\n".
							"\tSequence length = $lseq\n\tStructure length = $lstruc\n".
							"Check: Do these symbols, $str_f, match those in this file $file?\n";
				}
				($seqstring) or die "ERROR IN $FUNCTION: The sequence for header id $id is empty!\n".
									"Check: Do your sequences contain letters other than [ATGCU],\n?".
									"or is the format of your fasta file correct?\n";
				# check whether the format is known and if it fits to the structure given
				check_str_format($structure, $str_f, $FUNCTION);
				$fasta{$id}->{SEQUENCE} = $seqstring;
				$fasta{$id}->{STRUCTURE} = $structure;
				$seqstring 	= 	"";
				$structure 	= 	"";
			}
			$id = $1;
			push(@order, $id);
		} else {
			uc($line);
			if($line 		=~ 	/^\s*([AGTCU]+)\s*$/){
				$seqstring .= 	$1;
			}
			if($line 		=~ 	/^\s*([$str_f]+)\s*$/){
				$structure .= 	$1;
			} 
		}
	}
	
	if($id && $seqstring && $structure){
				my $lseq 	= 	length($seqstring);
				my $lstruc 	= 	length($structure);
				if($lseq   != 	$lstruc){
					die "ERROR IN $FUNCTION: The sequence length is not equal to the structure length in $file!\n".
							"\tSequence length = $lseq\n\tStructure length = $lstruc\n".
							"Check: Do these symbols, $str_f, match those in this file $file?\n";
				}
				($seqstring) or die "ERROR IN $FUNCTION: The sequence for header id $id is empty!\n".
									"Check: Do your sequences contain letters other than [ATGCU],\n?".
									"or is the format of your fasta file correct?\n";
				# check whether the format is known and if it fits to the structure given
				check_str_format($structure, $str_f, $FUNCTION);				
				$fasta{$id}->{SEQUENCE} = $seqstring;
				$fasta{$id}->{STRUCTURE} = $structure;
				
				$seqstring 	= 	"";
				$structure 	= 	"";
	}
	return (\%fasta, \@order);
}


##################################################################################
# Converts your input structure from the given format to the other given format,
# e.g. converts dot-bracket format to the constraint format and vice-versa.
#
# INPUT:
#		struct		The structure string
#		from_str	The structure format of the structure string
#		to_str		The structure format you want to convert the stucture into
#		check		(Optional) Check the given input formats
# OUTPUT:
#		The new structure sequence 
#
##################################################################################
sub convert_structure_format{
	my ($struct, $from_str, $to_str, $check) = @_;
	my $FUNCTION = "convert_structure_format in Structure.pm";
	
	die "INPUT ERROR in $FUNCTION: check your input variables!\n" 
		unless ($struct,$from_str, $to_str);
	
	if($check){
		$from_str = check_str_format($struct, $from_str, $FUNCTION);
		$to_str = check_str_format(0, $to_str, $FUNCTION);
	}
	
	# from dotbracket
	if($from_str eq $CONFIG{STR_FORMATS}->{DOTBRACKET}){
		# to constraint format
		if($to_str eq $CONFIG{STR_FORMATS}->{CONSTRAINT}){
			$struct =~ tr/\-\./\.x/d;
		# no other format implemented yet
		} else {
			die "ERROR in $FUNCTION: This structure format, $to_str, ".
				"has not been implemented yet!\n";			
		}
		
	# from constraint
	} elsif ($from_str eq $Config::CONFIG{STR_FORMATS}->{CONSTRAINT}){
		if($to_str eq $CONFIG{STR_FORMATS}->{DOTBRACKET}){
			$struct =~ tr/\.|<>x/\-\-\-\-\./d;
		} else {
			die "ERROR in $FUNCTION: This structure format, $to_str, ".
				"has not been implemented yet!\n";
		}
		
	} else {
		die "ERROR in $FUNCTION: This structure format, $from_str, ".
			"has not been implemented yet!\n";
	}
	
	return $struct;
}

##################################################################################
# Takes a dotplot and converts it into Fabrizio's graph format. The file contents 
# is returned as an array reference. The format is as follows:
# t # IDENTIFIER
# v VERTEX_NUMBER VERTEX VALUE
# e VERTEX1 VERTEX2 EDGE_SYMBOL EDGE_WEIGHT
# NOTE: The sequence order is given by an edge backbone and a special edge symbol.
# Also the numbering is from 0 to n-1.
#
# INPUT:
#		dotplot		(Compulsory) The dotplot file
#		graph_identifier	(Compulsory) The name of the graph with all information
#		edge_symbol	(Compulsory) The symbol for the edge weight
#		ts_seq		(Compulsory) The sequence of the target site, or region of interest
#		cutoff		(Optional) The cutoff value for the dotplot probabilities
# OUTPUT:
#		The graph as an array reference where each element is a row in the
#		output graph file for ONE graph
##################################################################################
sub convert_dotplot_to_fabrizio_graph{
	my ($dotplot, $graph_identifier, $edge_symbol, $ts_seq, $cutoff) = @_;
	my $FUNCTION = "convert_dotplot_to_fabrizio_graph in Structure.pm";
	
	$cutoff = 0 unless ($cutoff > 0);
	
	die "ERROR in $FUNCTION: The identifier parameter is compulsory!\n" unless ($graph_identifier);
	die "ERROR in $FUNCTION: You must give the symbol for the edge weight!\n" unless ($edge_symbol);
	die "ERROR in $FUNCTION: Symbol not allowed $edge_symbol for the edge!\n" if 
	($edge_symbol eq "#" || $edge_symbol eq "v" || $edge_symbol eq "e" || $edge_symbol eq ">");
	die "ERROR in $FUNCTION: You must give the target sequence!\n" unless ($ts_seq);
	
	my @graph; # contents of the graph format
	# read dotplot and extract base-pair probabilities and the sequence
	my $res_aref = parse_dotplot_return_bp_probs_and_sequence($dotplot);
	my $bp_href = $res_aref->[0];
	my $sequence = $res_aref->[1];
	
	#check for contents in the dotplot or if parsing produced a result
	die "ERROR in $FUNCTION: The sequence is empty in $dotplot, something went wrong\n" unless ($sequence);
	die "ERROR in $FUNCTION: There are no base pair probabilities in $dotplot\n" unless ($bp_href);
	
	# find ts position
	my $ts_index;
	$ts_seq  =~ tr/T/U/d;
	if (($ts_index = index($sequence, $ts_seq)) == -1){
		die "ERROR in $FUNCTION: The target sequence could not be found in the dotplot ".
			"sequence!\nTarget sequence: $ts_seq\nDotplot sequence: $sequence\n";
	}
	
	# add graph identifying line
	my $ts_len = length($ts_seq);
	push(@graph, "t # $graph_identifier TSpos=$ts_index-".($ts_index + $ts_len - 1));
	
	# add sequence information to graph
	my @seq_a = split('', $sequence);
	my $n = @seq_a;
	
	# add vertices
	for (my $i = 0 ; $i < $n ; $i++){
		push(@graph, "v $i $seq_a[$i]");
	}
	# add backbone edges, i.e. sequence order edges
	for (my $i = 0; $i < $n - 1 ; $i++){
		push(@graph, "e $i ".($i+1)." > 1");
	}
	
	# add base-pair probabilities as edges
	foreach (keys (%{$bp_href})){
		my ($i, $j) = split(":", $_);
		push (@graph, "e ".($i-1)." ".($j-1)." $edge_symbol $bp_href->{$_}") 
				if ($bp_href->{$_} >= $cutoff);
	}
	die "ERROR in $FUNCTION: The graph could not be built!\n" unless (@graph);
	return \@graph;
}

##################################################################################
# Takes a dotplot and adds it to an existing graph array reference. This means it
# checks whether it fits to the previous graph and then adds the edge weights of
# the new dotplot with the given edge symbol. The format is as follows:
# t # IDENTIFIER
# v VERTEX_NUMBER VERTEX VALUE
# e VERTEX1 VERTEX2 EDGE_SYMBOL EDGE_WEIGHT EDGE_SYMBOL EDGE_WEIGHT
# NOTE: The sequence order is given by an edge backbone and a special edge symbol.
# Also the numbering is from 0 to n-1.
#
# INPUT:
#		graph_aref	(Compulsory) The initial graph as an array reference
#		dotplot		(Compulsory) The dotplot file
#		graph_identifier	(Compulsory) The name of the graph with all information
#		edge_symbol	(Compulsory) The symbol for the edge weight
#		ts_seq		(Compulsory) The sequence of the target site, or region of interest
#		cutoff		(Optional) The cutoff value for the dotplot probabilities
# OUTPUT:
#		The updated graph as an array reference where each element is a row in the
#		output graph file for ONE graph
##################################################################################
sub add_dotplot_to_fabrizio_graph{
	my ($graph_aref, $dotplot, $graph_identifier, $edge_symbol, $cutoff) = @_;
	my $FUNCTION = "add_dotplot_to_fabrizio_graph in Structure.pm";
	
	$cutoff = 0 unless ($cutoff > 0);
	
	# check input parameters
	die "ERROR in $FUNCTION: The graph array reference is compulsory!\n" unless ($graph_aref);
	die "ERROR in $FUNCTION: The identifier parameter is compulsory!\n" unless ($graph_identifier);
	die "ERROR in $FUNCTION: You must give the symbol for the edge weight!\n" unless ($edge_symbol);
	die "ERROR in $FUNCTION: Symbol not allowed $edge_symbol for the edge!\n" if 
	($edge_symbol eq "#" || $edge_symbol eq "v" || $edge_symbol eq "e" || $edge_symbol eq ">"
	|| $edge_symbol eq "t");
	
	my @graph; # contents of the graph format
	# read dotplot and extract base-pair probabilities and the sequence
	my $res_aref = parse_dotplot_return_bp_probs_and_sequence($dotplot);
	my $bp_href = $res_aref->[0];
	my $sequence = $res_aref->[1];
	
	#check for contents in the dotplot or if parsing produced a result
	die "ERROR in $FUNCTION: The sequence is empty in $dotplot, something went wrong\n" unless ($sequence);
	die "ERROR in $FUNCTION: There are no base pair probabilities in $dotplot\n" unless ($bp_href);

	# read initial graph into a hash reference and check if we are
	# adding to the correct graph
	my @seq_a = split('', $sequence);
	my %graph;
	foreach (@{$graph_aref}){
		if ($_ =~ /^t #\s(.+)\sTSpos.*$/){
			die "ERROR in $FUNCTION: The initial graph array differs to the identifier ".
			"given!\nInitial identifier - $1\nGiven identifier - $graph_identifier\n"
			 	unless ($1 eq $graph_identifier);
			 push(@graph, $_);
		} elsif ($_ =~ /^v\s(\d+)\s([ATGCU])/){
			die "ERROR in $FUNCTION: The sequence in $dotplot does not match to the ".
				"original graph given!\n" unless ($seq_a[$1] eq $2);
			push(@graph, $_);
		} elsif ($_ =~ /^e\s\d+\s\d+\s\>\s1/){
			# this is the backbone and it remains the same
			push(@graph, $_);
		} elsif ($_ =~ /^e\s(\d+)\s(\d+)\s.+/){
			if(exists($bp_href->{"$1:$2"})){
				push(@graph, $_." $edge_symbol ".$bp_href->{"$1:$2"});
				delete($bp_href->{"$1:$2"});
			} else{
				# no base pair in new dotplot
				push(@graph, $_);
			}
		}
	}
	
	# push remaining base-pairs that did not exist in the previous graph
	foreach (keys (%{$bp_href})){
		my ($i, $j) = split(":", $_);
		push (@graph, "e $i $j $edge_symbol $bp_href->{$_}");
	}
	
	die "ERROR in $FUNCTION: The graph could not be built!\n" unless (@graph);
	return \@graph;
}

##################################################################################
# Calls Rfold from the paper Kiryu et al, Bioinformatics 2008. This program
# does a local folding, similar to the RNAlfold algorithm. So it calculates
# exact base-pair probabilities while only allowing a maximal distance between
# two base-pairs (parameter).
#
# INPUT:
#	fasta		The fasta file containing the input sequence(s)
#	dir			The output directory
#	L			max distance between two base-pairs
# OUTPUT:
#	none. (saves the results directly into the output dir given)
##################################################################################
sub call_rfold{
	my ($fasta, $dir, $L) = @_;
	
	my $FUNCTION = "call_rfold in Structure.pm";
	
	die "INPUT ERROR in $FUNCTION: no fasta file was given!\n" unless ($fasta);
	die "INPUT ERROR in $FUNCTION: the fasta file doesn't exist!$fasta\n" 
		unless (-e $fasta);
	die "INPUT ERROR in $FUNCTION: no output dir was given!" unless ($dir);
	die "INPUT ERROR in $FUNCTION: the output dir doesn't exist!$dir" unless (-e $dir);
	die "INPUT ERROR in $FUNCTION: no parameter L was given for max bp length!\n"
		unless ($L);
	
	# /home/maticzkd/cluster/src/rfold-0.1-1/src/run_rfold -max_pair_dist=50 -print_prob=true X51902.fasta
	my $params = " -max_pair_dist=$L -print_prob=true $fasta";
	
	$dir = "." unless ($dir);
	
	my $rfold_call = "cd $dir ; ".$CONFIG{RFOLD}.$params;
	
	system($rfold_call);
	1;		
}

##################################################################################
# Calls Rfold for ONE sequence only and saves the base-pair probability file
# into the given directory. Also it calculates the accessibilites and saves
# these into the given directory as well. It deletes the file with the ending
# rfold_out.txt as this is not needed.
#
# INPUT:
#	fasta		The fasta file containing only one sequence
#	dir			The directory in which to save the base-pair probs and acc files
#	L			The max distance between two base-pairs
#	prefix		The file name prefix of the two output files
#
# OUTPUT:
#	The file names with absolute path of (1) the base-pair probabilities and
# 	(2) the accessibilities in an array.
##################################################################################
sub rfold_1seq_and_modify_output{
	my ($fasta, $dir, $L, $prefix) = @_;
	
	my $FUNCTION = "Structure::rfold_1seq_and_modify_output";
	return 0 unless (-e $fasta);
	
	$dir = "." unless ($dir);
	
	## check that there is only one sequence in the fasta file
	my @fasta_res = StructureLibrary::Sequence::read_fasta_file($fasta);
	my $n = @{$fasta_res[1]};
	die "ERROR in $FUNCTION: The fasta file contains more than one ".
		"sequence $fasta\n" if ($n != 1);
	my $seq_len = length($fasta_res[0]->{$fasta_res[1]->[0]});
	
	## call rfold
	my $tmpdir = $dir."/tmp_rfold";
	system ("mkdir $tmpdir");
	$tmpdir = abs_path($tmpdir)."/";
	call_rfold($fasta, $tmpdir, $L);
	
	# read tmp directory
	my $rfold_prob_file = "";
	opendir(IMD, "$tmpdir") || die("ERROR: Cannot open directory $tmpdir");
		my @files= readdir(IMD);
		foreach my $file (@files){
			## delete rfold_out.txt as it is not needed
			if($file =~ /rfold_out.txt$/){
				system("rm ".$tmpdir.$file);
			} elsif ($file =~ /(rfold_prob.txt)$/){
				$rfold_prob_file = $tmpdir.$1;
			} else {
				print "WARNING: what is this file???? $file\n" unless(
					$file eq "." || $file eq ".." );
			}
		}
	closedir(IMD); 
	
	## die if no prob file could be found
	die "ERROR in $FUNCTION: there is not probability file\n" 
		unless ($rfold_prob_file);
		
	## rename the prob file and delete tmp_dir
	my $new_prob_name = $dir."/$prefix.prob";
	$new_prob_name = abs_path($new_prob_name);
	system("mv $rfold_prob_file $new_prob_name");
	$rfold_prob_file = $new_prob_name;
	die "ERROR in $FUNCTION: the rfold prob outfile doesn't exist!\n"
		unless (-e $rfold_prob_file);
	system("rm -R $tmpdir");
		
	## read probability file 
	my @bp_probs = read_rfold_prob_file_into_basepairs_hash($rfold_prob_file);
	die "ERROR in $FUNCTION: There are too many (or no) sequences in ".
		"$rfold_prob_file!\n" if (@bp_probs != 1);
	my $bp_href = shift(@bp_probs);
	
	## convert probabilities into accessibilities and create file
	my $acc_file = $dir."/$prefix.acc";
	$acc_file = abs_path($acc_file);
	convert_bp_to_accfile($acc_file, $bp_href, $seq_len);
	
	## return the file names with absolute path 
	return ($rfold_prob_file, $acc_file);
}


##################################################################################
# Converts the given base-pair probabilities into accessibilities for each sequence
# position. The results are saved to a file and the accessibilites are returned.
# The numbering is from 1..n.
#
# INPUT:
#	file		The file name for the accessibilites
#	bp_href		A hash reference to the base-pair probabilities (these need to be
#				exact without averaging)
#	seq_len		You have to give the sequence length, so that the gaps can be
#				filled in for bases that are not involved in any base-pair probs.
#	nofile		Set this flag if you do not want a file to be saved!
#	OUTPUT:
#	A hash reference containing the accessibilites for each position.
##################################################################################
sub convert_bp_to_accfile{
	my ($file, $bp_href, $seqlen, $nofile) = @_;
	my $FUNCTION = "Structure::convert_bp_to_accfile";
	
	## check input
	die "INPUT ERROR in $FUNCTION: no file was given!\n" unless ($file);
	die "INPUT ERROR in $FUNCTION: no base-pair probability hash was given!\n" 
		unless ($bp_href); 
	die "INPUT ERROR in $FUNCTION: the sequence length must be given!\n"
		unless ($seqlen);

	## calculate that pos i is paired from base-pair probs
	my $paired_prob_href = prob_base_i_is_paired($bp_href);
	
	my @acc_contents = ();
	my %acc = ();
	
	## iterate through each position, now from 0 to n-1
	for (my $pos = 1 ; $pos <= $seqlen ; $pos ++){
		my $unpaired_prob;
		if (defined($paired_prob_href->{$pos})){
			$unpaired_prob = 1 - $paired_prob_href->{$pos};
		} else {
			$unpaired_prob = 1;
		}
		push(@acc_contents, "$pos\t$unpaired_prob") unless ($nofile);
		$acc{$pos} = $unpaired_prob;
	}
	## write to file unless this is not wanted
	StructureLibrary::Tools::write_file_from_array($file, \@acc_contents) 
		unless ($nofile);
		
	## return accessibility values;
	return \%acc;
}

##################################################################################
# Here the output base-pair probability file from Rfold is read and saved into
# a base-pair hash, which is returned as a reference. Rfold is capable of folding
# a multi-sequence fasta file. The result is saved in one file only, in the same
# order as given in the fasta file. For this reason an array of base-pair prob.
# hashes are returned.
# NOTE: numbering from 1 to n.
# ADDITIONAL NOTE: Rfold output includes errors as some base-pairs are
# listed twice. ....find out which entry is incorrect
# 
# INPUT: 
#	rfold_prob_file	The results file from Rfold
# OUTPUT:
#	An array of hash references for the base-pair probabilities of each input
# 	sequence in the order as was given in the original fasta file.
##################################################################################
sub read_rfold_prob_file_into_basepairs_hash{
	my ($rfold_prob_file) = @_;
	
	my $FUNCTION = "Structure::read_rfold_prob_file_into_basepairs_hash";
	
	my @all_res = ();
	my %current_bp_hash;
	
	## read rfold_prob_file into a hash
	open(IN_HANDLE, "<$rfold_prob_file") || 
		die "ERROR in $FUNCTION:\n Couldn't open file $rfold_prob_file/n";
	
	while(my $line = <IN_HANDLE>){
		chomp($line);
		
		## if it is the start of a new sequence, then start a new hash
		## and add the old one to all results
		if($line =~ /#left_pos\s+right_pos\s+prob\s+type/){
			if(%current_bp_hash){
				push(@all_res, \%current_bp_hash);
			}
			my %new_entry = ();
			%current_bp_hash = %new_entry;
		}
		
		##if the line contains a base-pair probability add it to the hash
		elsif($line =~ /(\d+)\s(\d+)\s(\d\.?\d*\S+)\sP/){
			next if ($3 > 1); ## just to be sure, throw out the large 'probs'
			my $key = "";
			if($1 < $2){
				$key = "$1:$2";
			} else {
				$key = "$2:$1";
			}
			## The Rfold indices are in the file twice and
			## only the first set is correct. In the second 
			## set there are numbers above 1.
			unless(defined ($current_bp_hash{$key})){
				$current_bp_hash{$key} = $3;
			}
			
		}
		else {
			print "This line is being ignored: <START>$line<END>\n" 
				unless($line =~ /^\s*$/);
		}
	}
	## add final sequence
	if(%current_bp_hash){
		push(@all_res, \%current_bp_hash);
	}
	
	close(IN_HANDLE);
	
	return @all_res;
}

##################################################################################
# This generates a dotplot that averages the base-pair probabilities of all
# repeat positions in a CRISPR array. This average dotplot gives a good overview
# of which base-pairs are most likely across the entire array. From this you
# can see the structure, which is most likely to be functional. Of course under
# the assumption that the functional struture is most prominent when folding
# the entire array.
# The top right side gives you the average base-pair probabilities with the
# "best structure" highlighted in red. The bottom left half shows you the
# mfe structure of the direct repeat folded on its own.
# 
# INPUT: 
#	bp_prob_href	This is the base-pair probability hash reference from the 
#					reference dotplot (using method parse_dotplot*)
#	array_seq		The sequence of the entire CRISPR array
#	dr_seq			The sequence of the direct repeat in the array
#	best_structure	(Optional) This is the "best" structure, according to 
#					some evaluation criteria. E.g. the maximum average base-pair
#					accuracy across all repeat positions
#	dr_mfe			This is the MFE structure of the direct repeat
#	str_format		The format of the best_structure parameter see 
#					$CONFIG{STR_FORMATS}
#	avedot_file		This is the name of the output dotplot file
#	
# OUTPUT:
#	none
##################################################################################
sub average_crispr_dotplot{
	my ($bp_prob_href, $array_seq, $dr_seq, $best_structure, $dr_mfe,
		$str_format, $avedot_file) = @_;
	my $FUNCTION = "Structure::ave_dotplot_of_CRIPSR_repeats";
	
	my $dr_len = length($dr_seq);
	my @positions = StructureLibrary::Sequence::find_motif($array_seq, $dr_seq);
	(@positions) or die "ERROR in $FUNCTION: the repeat sequence could not be ".
		"found in the array sequence!\n"; 
	my $num_p = @positions;
	
	die "ERROR in $FUNCTION: the repeat sequence does not occur in the array!\n"
		unless ($num_p);
		
	## fill base-pair probabilities with 0
	my %ave_bp_prob = ();
	for (my $i = 1 ; $i <= $dr_len ; $i ++){
		## no loops the size of 0,1,2!
		for (my $j = $i+3 ; $j <= $dr_len ; $j++){
			$ave_bp_prob{"$i:$j"} = 0;
		}
	}
	
	
	## calculate average base-pair probabilites
	
	foreach my $p (@positions){
		for (my $i = 1 ; $i <= $dr_len ; $i ++){
			## no loops the size of 0,1,2!
			for (my $j = $i+3 ; $j <= $dr_len ; $j++){
				my $key = ($p + $i).":".($p+$j);
				if(exists($bp_prob_href->{$key})){
					$ave_bp_prob{"$i:$j"} += $bp_prob_href->{$key}/$num_p;
				} else {
					$ave_bp_prob{"$i:$j"} += 0;
				}
			}
		}
	}
	
	## get mfe structure co-ordinates
	die "ERROR in $FUNCTION: no best structure given!\n" unless ($dr_mfe);
	my $mfe_bp_href = structure2basepairs(
			$dr_mfe, 1, $str_format);
	
	
	## generate template dotplot
	system ("echo \"$dr_seq\" | $CONFIG{RNAFOLD} $CONFIG{RNAFOLDPARAMS}");
	system ("rm rna.ps");
	system ("mv dot.ps $avedot_file");
	
	## get the contents of the templat that is needed
	my @dot_contents = StructureLibrary::Tools::read_file_return_array($avedot_file);
	my @ave_dot_contents = ();
	
	foreach my $line (@dot_contents){
		chomp($line);
		push(@ave_dot_contents, $line);
		if($line =~ /\%data\sstarts\shere/){
			last;
		}
	}
	
	## fill ubox dots
	foreach (keys %ave_bp_prob){
		my ($i,$j) = split(":", $_);
		push(@ave_dot_contents, "$i $j ".$ave_bp_prob{$_}." ubox");
	}
	
#	for (my $i = 1 ; $i <= $dr_len ; $i ++){
#		## no loops the size of 0,1,2!
#		for (my $j = $i+3 ; $j <= $dr_len ; $j++){
#			push(@ave_dot_contents, "$i $j ".$ave_bp_prob{"$i:$j"}." ubox");
#		}
#	}
	
	## fill lbox dots (mfe structure)
	foreach (keys (%{$mfe_bp_href})){
		my ($i,$j) = split(":", $_);
		push(@ave_dot_contents, "$i $j 1 lbox");
	}
	
	## end ps file
	push(@ave_dot_contents, "showpage");
	push(@ave_dot_contents, "end");
	push(@ave_dot_contents, "%%EOF");
	
	## write file
	StructureLibrary::Tools::write_file_from_array($avedot_file, \@ave_dot_contents);
	
	## highlight best structure in dotplot
	if($best_structure){
		my $ofile = $avedot_file."_marked.pdf";
		if($avedot_file =~ /(.+)\.ps/){
			$ofile = $1."_marked.pdf";
		} 
		my $call = $CONFIG{DOTVIEW}." -dot $avedot_file -o $ofile -str-db \"".
			"red, 1, $best_structure\" -pdf";
		system($call);
	}
}

##################################################################################
# Taking a base-pair hash reference, this method calculates the maximum span
# between the base-pairs. The max span of the given structure is the max. 
# distance between a base-pair (i,j), which is exactly |j-i| +1, meaning
# the distance between i and j on the sequence including positions i and j.
#
# INPUT: 
#	bp_href			(Compulsory)
#					This is the hash reference including all base-pairs for the 
#					given structure. 
#	
# OUTPUT:
#	The max. span of the given structure.
##################################################################################
sub get_maxspan_of_structure{
	my ($bp_href) = @_;
	my $FUNCTION = "Structure::get_maxspan_of_structure";
	
	## check input
	die "INPUT ERROR in $FUNCTION: bp_href is missing!\n" unless $bp_href;
	
	## get maxspan
	my $maxspan = 0;
	foreach my $bps (keys (%{$bp_href})){
		my ($i,$j) = split(":", $bps);
		my $span = abs($j - $i) + 1; # num nucleotides including i and j
		$maxspan = $span if ($span > $maxspan);
	}
	die "ERROR in $FUNCTION: the maxspan is <4, this must be an error!\n" 
		if ($maxspan < 4);
	return $maxspan;
}


1;