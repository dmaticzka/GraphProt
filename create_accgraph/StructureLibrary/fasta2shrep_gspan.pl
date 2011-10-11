#!/usr/bin/perl
#use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use StructureLibrary::RNAtoolsConfig;
use StructureLibrary::Tools;
use StructureLibrary::Sequence;
use StructureLibrary::Structure;
use Cwd qw(getcwd abs_path);


use vars qw(%CONFIG);
*CONFIG = \%StructureLibrary::RNAtoolsConfig::CONFIG;


=head1 NAME


fasta2shrep_gspan.pl -fas mysequences.fasta -wins "50,100,150,200" -shift 5 -M 8

=head1 SYNOPSIS

Options:

		HELP
        -help   brief help message
        -man    full documentation
        
        COMPULSORY
        -fasta	<STRING> e.g. "sequence.fasta"
        		All sequences in fasta format.
        		
       	OPTIONS
        -wins	[INTEGER] e.g. "50,100,200"
        		A list of window sizes to use.
        		If none are given, then the entire sequence is taken with no windows.
        -shift	<INTEGER> e.g. 5
        		The number of nucleotides the window should be shifted 
        		by for each calculation (RNAshapes)
        		If none is given, the default is a shift of 1.
        -cue	Crop unpaired ends.
        		If you give this flag, then the unpaired ends of each
        		single structure are ignored. E.g. the structure
        		...(((...))).. becomes just (((...)))
        -e		<FLOAT> e.g. 5.0
        		Energy range in kcal/mol (RNAshapes)
        		Use only one of -e and -c!
        -c 		<INTEGER> e.g. 10
        		Relative energy range, i.e. percentage (%) of MFE energy (RNAshapes)
        		Use only one of -e and -c!
        -t		<INTEGER> [1-5] e.g. 3
        		The shape type (RNAshapes). Default is 3.
        -M		<INTEGER> e.g. 10
        		Max number of shreps that should be taken per window.
        -u 		Ignore unstable structures (RNAshapes). 
        		This option filters out closed structures with positive free energy. 
        -r		Calculate structure probabilities for shreps (RNAshapes)
        -tmp	<STRING> e.g. "/scratch/1/sita/tmp"
        		A directory for writing temporary files
        -o		<STRING> e.g. "ProjectX/MySequences/GSPAN/" 
        		Output directory for gspan files containing graphs.
        -sge	Use SGE cluster for each sequence separately
        		
        		
        DEFAULT VALUES
        -wins	""
        -shift	5
        -c		10
        -t		3
        -M		10
        -tmp	"/var/tmp/fasta2shrep"
        -o		"CURRENT_DIR/GSPAN/"


=head1 DESCRIPTION

=cut


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################

# command line options
my ($i_help, $i_man, $i_fas, $i_wins, $i_shift, $i_crop_unpaired_ends,$i_r, 
$i_e, $i_c, $i_t, $i_u, $i_M,  $i_o, $i_sge, $i_jobid, $i_tmp);


my $options = GetOptions ("help"	=> \$i_help,
                         "man"  	=> \$i_man,
                         "fasta=s"	=> \$i_fas,
                         "wins=s"	=> \$i_wins,
                         "shift=i"	=> \$i_shift,
                         "cue"		=> \$i_crop_unpaired_ends,		
                         "r"		=> \$i_r,
                         "e=f"		=> \$i_e,
                         "c=i"		=> \$i_c,
                         "t=i"		=> \$i_t,
                         "u"		=> \$i_u,
                         "M=i"		=> \$i_M,
                         "tmp=s"	=> \$i_tmp,
                         "o=s"		=> \$i_o,
                         "sge"		=> \$i_sge,
                         "jobid=s"	=> \$i_jobid
                         );

                      
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;
pod2usage(-exitstatus => 0, -verbose => 2) if $i_man;
($options) or pod2usage(2);

# check compulsory options
($i_fas) or pod2usage("Error: the option -fas is compulsory!\n");
(-e $i_fas) or pod2usage("Error: no such file - $i_fas!\n");
$i_fas = abs_path($i_fas);

# check other options and set default values
pod2usage("Error: either -e OR -c can be given, but NOT both!\n") if($i_e && $i_c);
($i_e or $i_c) or $i_c = 10; # set -c 10, if neither -e or -c are given
($i_t) or $i_t = 3; # default abstraction type is 3
($i_M) or $i_M = 0; # max number of shreps is 10
my $CURRDIR = getcwd;
if($i_tmp){
	(-e $i_tmp) or system("mkdir -p $i_tmp");
} else {
	## default tmp is /var/tmp, usually not NFS mounted!
	system("mkdir -p /var/tmp/fasta2shrep");
	$i_tmp = "/var/tmp/fasta2shrep";
}
if($i_o){
	(-e $i_o) or system("mkdir -p $i_o");
} else {
	system("mkdir -p GSPAN");
	$i_o = $CURRDIR."/GSPAN/";
}

if($i_wins){
	($i_shift) or pod2usage("Error: the option -shift is compulsory if -wins is given!\n");
} else {
	pod2usage("Error: the option -shift must be given with -wins!\n") if ($i_shift);
}

###############################################################################
# GLOBAL VARIABLES
###############################################################################
my $SCRIPTDIR = $CONFIG{LIB_DIR};
my $SCRIPTNAME = "$SCRIPTDIR/fasta2shrep_gspan.pl";
my $CLUSTERSUBMIT = "$SCRIPTDIR/generic_submit_to_cluster.sh";
my @WINDOWS = ();
@WINDOWS = split(",", $i_wins) if ($i_wins);
my $globalFolding;
$globalFolding = 1 unless @WINDOWS;
my $CURRUSER = getlogin;


###############################################################################
# EXECUTION CODE
###############################################################################

# read fasta file into hash
my $in_fasta = read_fasta_file($i_fas);
my $sequences_href = $in_fasta->[0];
my $seq_order_aref = $in_fasta->[1];
my $seq_heads_href = $in_fasta->[2];

# call script again on the sge cluster in a batch job
if($i_sge){
	my $params = "-fasta $i_fas ";
	$params .="-wins $i_wins " if ($i_wins);
	$params .= "-shift $i_shift " if ($i_shift);
	$params .= "-e $i_e " if ($i_e); 
	$params .= "-c $i_c " if ($i_c);
	$params .= "-t $i_t " if ($i_t);
	$params .= "-u " if ($i_u);
	$params .= "-r " if ($i_r);
	$params .= "-M $i_M " if ($i_M);
	$params .= "-o $i_o " if ($i_o);
	$params .= "-tmp $i_tmp " if ($i_tmp);
	
	my $ssh=1;
	my $mRNA_num = @{$seq_order_aref};
		
	if ($ssh){
		system("ssh $CURRUSER\@biui.informatik.uni-freiburg.de ".
			"'export SGE_ROOT=/opt/sge-6.0/; cd $CURRDIR; ".
			"/opt/sge-6.0/bin/lx24-amd64/qsub -t 1-$mRNA_num ".
			"$CLUSTERSUBMIT $CURRDIR \"$SCRIPTNAME\" \"$params\" ' ");
	}else{
		print "params: $params\n";
		system("$CLUSTERSUBMIT $CURRDIR \"$SCRIPTNAME\" \"$params\" ");
	}	
		
} else {
	
	my @used_seq_ids;
	
	if ($i_jobid){
		## just process the one sequence as given by the jobid number
		push(@used_seq_ids,$seq_order_aref->[$i_jobid -1]);
	} else{
		## process all sequences at once
		@used_seq_ids = sort keys %{$sequences_href};
	}
		
	# for each sequence in the fasta file (any order)
	foreach my $seq_id (@used_seq_ids){

		my $seq_fasta = generate_single_fasta_from_sequence_X($seq_id,$sequences_href->{$seq_id});
		my $seq_header = $seq_heads_href->{$seq_id} if (exists $seq_heads_href->{$seq_id});
			
		if ($globalFolding){
			@WINDOWS = ();
			push(@WINDOWS,length($sequences_href->{$seq_id}));
		}

		my $gspanfile = $i_o."/$seq_id.gspan";
		open(OUT,"| bzip2 -f > $gspanfile.bz2");
			
		my $graph_header = getGraphHeader($seq_id,$seq_header,\@WINDOWS,$i_shift,$i_e, $i_c, $i_t, $i_u, $i_r, $i_M);
		print OUT $graph_header;
		
		my $gi = 0; # current graph index
		
		#for each window size in list
		foreach my $win_size (@WINDOWS){
				
			print "Next: $seq_id\t winsize:$win_size\n";
				
			my $rnashapes_outfile = $i_tmp."/$seq_id.W${win_size}.RNAshapes.out";
				
			# call RNAshapes and write to file
			call_RNAshapes($rnashapes_outfile, $seq_fasta, 
							$CONFIG{RNASHAPES}, $win_size, $i_shift,
							$i_e, $i_c, $i_t, $i_u, $i_r);
		
			# read RNAshapes output and write subgraph to gspan file
			$gi = convert_RNAshapes_output($rnashapes_outfile, $gi,$i_M,\*OUT,$graph_header);
				
			# remove RNAshapes output file as we don't need it anymore
			system("rm $rnashapes_outfile");
		}
		close(OUT);
	}
}

###############################################################################
# METHODS
###############################################################################

sub generate_single_fasta_from_sequence_X{
	my ($seq_id, $seq) = @_;
 
	$seq = uc($seq);
	$seq =~ tr/T/U/d;
	$seq =~ s/[^AUCGN]/N/g;
	
	my $outfas = $i_tmp."/$seq_id.fasta";
	
	open(FAS,">$outfas") or die "Cannpt open file $outfas! Exit...\n\n";
	print FAS ">$seq_id\n$seq";
	close(FAS);
	
	return $outfas;
}

sub call_RNAshapes{
	my ($rnaoutfile, $seq_fasta, $rnashapes_location, $win_size, $shift, $e, $c, 
		$t, $u, $r) = @_;
	my $FUNCTION = "call_RNAshapes in fasta2shrep_gspan.pl";
	
	($rnaoutfile) or die("INPUT ERROR in $FUNCTION: the output file name for ".
						"RNAshape is compulsory!\n");
	($seq_fasta) or die("INPUT ERROR in $FUNCTION: the fasta file is compulsory!\n");
	($rnashapes_location) or die ("INPUT ERROR in $FUNCTION: the RNAshapes location".
									" is compulsory!\n");
	
	my $call = $rnashapes_location." -o 1 "; # the output format is of type 1
	$call .= "-w $win_size "; 
	$call .= "-W $shift " if ($i_shift);
	die("ERROR in $FUNCTION: Give only one of the options -c or -e (RNAshapes)!\n")
		if ($e && $c);
	$call .= "-e $e " if ($e); 
	$call .= "-c $c " if ($c);
	$call .= "-t $t " if ($t);
	$call .= "-u " if ($u);
	$call .= "-r " if ($r);
	$call .= " < $seq_fasta > $rnaoutfile";
	
	system ($call);
	
	(-e $rnaoutfile) or die("ERROR in $FUNCTION: The following call ".
							"could not be carried out!\n$call\n");
	
	1;
}

sub convert_RNAshapes_output{
	my ($rnashapesoutput, $curr_gi, $maxShreps, $graph_file_hdl, $graphHead) = @_;
	
	open(IN,"$rnashapesoutput");

	## seqId, not used
	my $line = <IN>;
	my @win_shreps;
	my $win_start;
	my $win_end;
	my $win_seq;
	my $win_shrep_count = 0;

	while ($line = <IN>){
		
		if ($line =~ /^(\d+)\s+(\d+)$/){
			## line: "<start>    <end>"
			if (@win_shreps > 0){
				#print " Found shreps=".@win_shreps."\n";
				$curr_gi = convertShapeWindow(\@win_shreps,$win_seq,$win_start,$win_end,$curr_gi,$graph_file_hdl,$graphHead)
			}
			
			## set new window params
			@win_shreps = ();
			$win_start = $1;
			$win_end = $2;
			$win_shrep_count = 0;

			#print "Window ".($win_end-$win_start+1)." ($win_start -$win_end):";
			
		} elsif ($line =~ /^(\S+)$/){
			## line: "CUUAUGAGUAAGGAAAAUAACGAUUCGGGGUGACGCCCGAAUCCUCACUG"
			$win_seq = $1;
	
		} elsif ($line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+(\S+)$/){
			## line:"...((((..(((....)))))))...........(((((......)))))  (-10.10)  [[]][]"
			## take only $maxShreps shreps per window if set
			next if ($maxShreps && $win_shrep_count>=$maxShreps);
			push(@win_shreps,[$1,"ENERGY",$2,"SHAPE",$3]);
			$win_shrep_count++;
			
		} elsif ($line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+\((\S+)\)\s+(\S+)$/){
			## line:"((((..((((...((.((.((.....)).)).))...))))..))))...  (-10.60)  (0.7795360)  [[[[[]]]]]"
			## take only $maxShreps shreps per window if set
			next if ($maxShreps && $win_shrep_count>=$maxShreps);
			push(@win_shreps,[$1,"ENERGY",$2,"PROB",$3,"SHAPE",$4]);
			$win_shrep_count++;

		} else {
			die "Unexpected shape output format!\nline=$line\n\nExit...\n\n";
		}	
	}
	
	## convert last window
	if (@win_shreps > 0){
		#print " Found shreps = ".@win_shreps."\n";
		$curr_gi = convertShapeWindow(\@win_shreps,$win_seq,$win_start,$win_end,$curr_gi,$graph_file_hdl,$graphHead)
	}
	
	close(IN);
	
	return $curr_gi + 1; # return the gi (graph index) for the next subgraph
}

sub convertShapeWindow{
	my ($win_shreps_aref,$win_seq, $win_start,$win_end,$curr_gi,$graph_file_hdl,$graphHead) = @_;
	 
	my $win_size = $win_end - $win_start + 1;
	my $win_center = $win_start+(($win_size+1)/2);

	my $winHead = getWindowHeader($graphHead,$win_size,$win_start,$win_end, $win_center);
	print $graph_file_hdl $winHead;
	
	## generate for each shrep a connected component in gspan file
	foreach my $shrep (@{$win_shreps_aref}){
		
		# cut off unpaired ends, if option is given
		my $shrep_seq = $shrep->[0];
		my $crop_index_left = 0;
		my $crop_index_right = length($shrep_seq) -1;
		
		if($i_crop_unpaired_ends){
			
			$crop_index_left = index($shrep_seq, "("); # find 1st occ of "("
			$crop_index_right = rindex($shrep_seq, ")"); # find last occ of ")"	
			# if the complete window is unpaired, then don't crop
			if($crop_index_left == -1){ 
				$crop_index_left = 0;
				$crop_index_right = length($shrep_seq) - 1;
			}		
		}
		
		my $backboneGraph_ref = getBackboneGraph($win_seq,$curr_gi,$win_start, 
								$crop_index_left, $crop_index_right);
		my $structGraph_ref = getStructGraph($shrep,$curr_gi);

		print $graph_file_hdl getShapeHeader($winHead,$shrep);
		print $graph_file_hdl join("\n",@{$backboneGraph_ref})."\n";
		print $graph_file_hdl join("\n",@{$structGraph_ref})."\n";
		
		$curr_gi = $curr_gi + $win_size;
	} 
	
	return $curr_gi;
}

sub getBackboneGraph{
	my ($win_seq,$curr_gi,$win_start, $curr_crop_i_left, $curr_crop_i_right) = @_;

	my @seq = split("",$win_seq);
	
	my @vert	= ();
	my @edg		= ();

	my $curr_abs_pos = $win_start+$curr_crop_i_left;
	$curr_gi += $curr_crop_i_left;
	
	push(@vert,"v $curr_gi $seq[$curr_crop_i_left] $curr_abs_pos");

	foreach my $idx (($curr_crop_i_left+1)..$curr_crop_i_right){
		$curr_abs_pos++;
		$curr_gi++;
		push(@vert,"v $curr_gi $seq[$idx] $curr_abs_pos");
		push(@edg,"e ".($curr_gi-1)." ".$curr_gi." > 1");
	}
	
	my @ret = (@vert,@edg);
	
	return \@ret;
}

sub getStructGraph{
	my ($curr_shrep,$curr_gi) = @_;
	
	my @struct = split("",$curr_shrep->[0]);
	
	my @edg		= ();
	my @pairs	= ();

	foreach my $idx (0..@struct-1){
		
		push(@pairs,$idx) if ($struct[$idx] eq "(");
		
		if ($struct[$idx] eq ")"){
			my $start = pop(@pairs);
			push (@edg,"e ".($curr_gi+$start)." ".($curr_gi+$idx)." s ");
		}
	}
	
	return \@edg;
}
	
sub getGraphHeader{
	my ($seq_id,$seq_head,$wins_aref,$h_shift,$h_e, $h_c, $h_t, $h_u, $h_r, $h_M) = @_;
	
	my $ret;
	
	$ret  = "t # SEQID $seq_id ";
	$ret .= "$seq_head " if ($seq_head); 
	$ret .= "MAXWINSHREPS $h_M ";
	$ret .= "RNASHAPES -w ".(join(",",@{$wins_aref}))." ";
	$ret .= "-W $h_shift " if ($h_shift);
	$ret .= "-e $h_e " if ($h_e); 
	$ret .= "-c $h_c " if ($h_c);
	$ret .= "-t $h_t " if ($h_t);
	$ret .= "-u 1 " if ($h_u);
	$ret .= "-r 1 " if ($h_r);
	$ret .= "\n";
	
	return $ret;
}

sub getWindowHeader{
	my ($graphHead,$win_size,$win_start,$win_end, $win_center) = @_;
	
	chomp $graphHead;
	
	$graphHead =~ s/t # //;
		
	my $ret = "w # $graphHead "; 
	$ret .= "WSIZE $win_size ";
	$ret .= "WSTART $win_start ";
	$ret .= "WEND $win_end ";
	$ret .= "WCENT $win_center";
	$ret .= "\n";
	
	return $ret;
}

sub getShapeHeader{
	my ($winHead,$shrep) = @_;
	
	chomp $winHead;
	
	$winHead =~ s/w # //;
		
	my $ret = "s # $winHead ";
	
	my @info = @{$shrep};
	$ret .= join(" ",@info[1..$#info]);	
	$ret .= "\n";
	
	return $ret;
}

1;

#############################################################################
# Programming description
#############################################################################
#	Substructure graphs for machine learning with Fabrizio
#	-------------------------------------------------------
#	
#	(1) Parameters (RNAshapes parameter):
#		- Window sizes [] (-w)
#		- window shift size (-W)
#	- calculate structure probabilities (-r)
#		- energy range kcal/mol (-e) OR energy relative percentage (%) to MFE (-c)
#		- shape type 1 to 5 (-t)
#		- ignore unstable substructures (-u)
#		- max shreps
#	
#	
#	(2) For each sequence, generate one graph/file that consists of all windows. The general format for one graph is as follows:
#	
#	t # seq_id parameters
#	v graph_index nt_type window_size window_centre abs_seq_pos
#	...
#	e m m+1 > 1 (backbone)
#	...
#	e base_i_graph_index base_j_graph_index s shrep_e shrep_p ...
#	
#	For each window (subgraph) we create a subgraph (of the subgraph) for each substructure.
#	We have a running index (gi) for each subgraph. All vertex and edge indices of the subgraph add
#	the running graph index to the actual window position. For example 
#	
#	
#	Sequence: AAACC CUUUG GG
#		  01234 56789 01
#	
#	Window=10 substructure1 = (((...))). centre 5.5
#	
#	v 0 A 10 5.5 0
#	v 1 A 10 5.5 1
#	v 2 A 10 5.5 2
#	v 3 C 10 5.5 3
#	v 4 C 10 5.5 4
#	v 5 C 10 5.5 5
#	v 6 U 10 5.5 6
#	v 7 U 10 5.5 7
#	v 8 U 10 5.5 8
#	v 9 G 10 5.5 9
#	e 0 1 > 1
#	e 1 2 > 1
#	e 2 3 > 1
#	e 3 4 > 1
#	e 4 5 > 1
#	e 5 6 > 1
#	e 6 7 > 1
#	e 7 8 > 1
#	e 8 9 > 1
#	e 0 8 s -15.0 0.1223
#	e 1 7 s -15.0 0.1223
#	e 2 6 s -15.0 0.1223
#	
#	gi = 9+1 = 10
#	
#	Window = 10 substructure2 = .(((...))) centre 6.5
#	v 10 A 10 6.5 2
#	v 11 C 10 6.5 3
#	v 12 C 10 6.5 4
#	v 13 C 10 6.5 5
#	v 14 U 10 6.5 6
#	v 15 U 10 6.5 7
#	v 16 U 10 6.5 8
#	v 17 G 10 6.5 9
#	v 18 G 10 6.5 10
#	v 19 G 10 6.5 11
#	e 10 11 > 1
#	e 11 12 > 1
#	e 12 13 > 1
#	e 13 14 > 1
#	e 14 15 > 1
#	e 15 16 > 1
#	e 16 17 > 1
#	e 17 18 > 1
#	e 18 19 > 1
#	e 11 19 s -17.0 0.156
#	e 12 18 s -17.0 0.156
#	e 13 17 s -17.0 0.156
#	
#	gi = 19+1 = 20
#	
#	
#	
#	Write one perl script to create graphs for a set of sequences, fasta2shrep_gspan.pl.
#	
#	INPUT: 
#		-f fasta file with all sequences to compute
#		parameters as above
#	
#	OUTPUT:
#		one file per sequence that contains graph called seq_id.gspan
#	
#	(1) for each window size call RNAshapes and write to a tmp file
#	(2) parse result of RNAshapes (catch RNAshapes error - sequence too long?) check for max shreps.
#	(3) convert RNAshapes result to subgraph -> write to file (readpipe) look at efficiency and errors
#	(4) repeat (1) to (3) for each sequence