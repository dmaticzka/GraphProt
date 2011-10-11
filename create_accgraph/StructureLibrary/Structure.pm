package StructureLibrary::Structure;
use Cwd 'abs_path';
use strict;
use StructureLibrary::RNAtoolsConfig; # Config.pm must be in the same directory

# global variables
# This declares the named variables as package globals in the current package.
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %CONFIG);
*CONFIG = \%RNAtoolsConfig::CONFIG; # Global Configuration now aliased to global %CONFIG


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
#		seq			(Compulsory) The original/reference sequence. 
#					The structure is checked to see if it feasible given
#					this sequence (i.e. can the base-pairs form and are the
#					indices within the correct range?). Also the indices of
#					the PU values are checked to be within the sequence boundaries.
#		seq_len		(Optional) Give the sequence length if it has been calculated
#					previously.
#		normalise	(Optional) Normalise the accuracy score according to the number of
#					probabilities in the sum
#		silent		(Optional) If you are very sure of what you are doing you can suppress all
#					outputs, such as warnings.
# Output:
#		The accuracy score.
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
	die("INPUT ERROR in $FUNCTION:\nYou must give the reference sequence!\n") unless ($seq);
	$seq_len = length($seq) unless ($seq_len);
	
	# accuracy score
	my $acc = 0;
	my $normalise_sum = 0;
	my $bpacc = 0;
	my $puacc = 0;
	
	# get base-pairs from structure	counting 1 to n
	$str_f = $CONFIG{STR_FORMATS}->{DOTBRACKT} unless ($str_f);
	my $bp_href = structure2basepairs($struct, $start, $str_f);
	
	# check structure
	check_bps_wrt_seq($seq, $bp_href);
	
	
	# add the probability of each base-pair in the structure and times it by 2
	foreach my $bp (keys (%{$bp_href})){
		if(exists($bp_prob_href->{$bp})){
			$bpacc += $bp_prob_href->{$bp} * 2;
			$normalise_sum += 2;
		}
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
		$pu_aref = parse_PUfile_get1u($pufile, 1, $num_header);
				
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
		++ $normalise_sum;
	}
	
	die ("ERROR in $FUNCTION: No structure probabilities have been added to the accuracy score!\n") unless ($normalise_sum);
	
	$acc = $bpacc + $puacc;
	$acc = $acc/$normalise_sum if ($normalise);
#	print STDERR "base pair contribution: $bpacc\nunpaired contribution: $puacc\nfinal acc: $acc\n";
	
	return $acc;
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
#		bps			A hash reference containing base-pairs as keys in the format "i:j"
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
				if($base_i eq "C" && !$base_j eq "G");
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
	return \@pu;	
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
#		i<j and empty values, so you can check whether a base-pair exists in the 
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
				$bp{(pop(@open)).":".($i+$start)} = 0;
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
#		ENSEMBLE_E	The ensemble energy of the sequence as written to STDOUT by RNAfold
#		CONSTRAINT_E The ensemble energy of all structures that fit to the constraint
#		MFE_E		The energy of the MFE structure 
#		ED			The energy difference between all structures and constraint structures
#		CONSTRAINT_P The probability of the constraint structures
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
	
	# add constraint ensemble energy to the results hash
	$res_all_href->{CONSTRAINT_E} = $res_constraint_href->{ENSEMBLE_E};
	
	# calculate ED
	my $ed = $res_all_href->{ENSEMBLE_E} - $res_all_href->{CONSTRAINT_E};
	$res_all_href->{ED} = $ed;
	
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
		$tmpdir = "tmp_rnafold/";
		
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
		my $bp_file = "";
		my $ss_file = "";
	
		# read tmp directory
		my $dp_file;
		opendir(IMD, ".") || die("ERROR: Cannot open directory $tmpdir");
			my @files= readdir(IMD);
			foreach my $file (@files){
				if($file =~ /.*\/tmp_rnafold\/(\S+_dp\.ps)$/){
					$dp_file = $1;
				} elsif ($file =~ /.*\/tmp_rnafold\/(\S+_ss\.ps)$/){
					$ss_file = $1;
				}
			}
		closedir(IMD); 
		
		# move result files and delete tmpdir
		die ("ERROR in $FUNCTION: could not move $dp_file!\n") unless(system("mv $dp_file ../."));
		die ("ERROR in $FUNCTION: could not move $ss_file!\n") unless(system("mv $ss_file ../."));
		die ("ERROR in $FUNCTION: could not delete $tmpdir!\n") unless(system("cd ../ ; rm -R $tmpdir"));
		
		# save files in hash
		$results{DP} = $dir.$dp_file;
		$results{SS} = $dir.$ss_file;
	}
	
	# get energies
	my ($e_all, $e_mfe) = "";
	foreach (@RNAfold_results) {
		chomp($_);
	    # get line with ensemble-energy 
	    if ($_ =~ /.+\[\s?(\-\d+\.\d+)\]$/) {
	       	 	$e_all = $1;
	    } elsif ($_ =~ /.+\(\s?(\-\d+\.\d+)\)$/){
	    		$e_mfe = $1;
	    }
	}
	
	# If the call above did not include the energies, then something went wrong!
	if(!$e_all || !$e_mfe){
		my $RNAout = "";
		foreach (@RNAfold_results){
			$RNAout .= $_."\n";
		}
		die "ERROR in $FUNCTION: see output below!\n $RNAout\n";
	}
	
	# add energies to hash and return it
	$results{ENSEMBLE_E} = $e_all;
	$results{MFE_E} = $e_mfe;
	return \%results;
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
			if($line =~ /^(\d+)\s(\d+)\s(\d\.\d+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			# modified dotplots
			} elsif($line =~ /rgb.*\s(\d+)\s(\d+)\s(\d\.\d+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			}
		}
		
	}
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
			if($line =~ /^(\d+)\s(\d+)\s(\d\.\d+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			# modified dotplots
			} elsif($line =~ /rgb.*\s(\d+)\s(\d+)\s(\d\.\d+)\subox/){
				if($1 < $2){
					$bp{"$1:$2"} = $3*$3 if ($3*$3 >0);
				# this case I expect never to occur, but it is handled just in case
				} else {
					$bp{"$2:$1"} = $3*$3 if ($3*$3 >0);
				}
			}
		}
	}
	
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
	
	my @result = (\%bp, $seqline);
	return \@result;
}



##################################################################################
# Reads the base-pair probabilites as returned by "parse_dotplot_return_bp_probs"
# and adds the probabilites for each base to a single sum. This sum is the
# probability for the respective base to be PAIRED. A hash is returned,
# where the key is the index for the base and the value is the sum of all 
# probabilites of base-pairs in which that base is involved in, i.e.
# the probability for the base to be paired.
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
	my @return = (\%fasta, \@order);
	return \@return;
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




1;