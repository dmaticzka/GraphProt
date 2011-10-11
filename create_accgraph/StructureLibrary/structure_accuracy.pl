#!/usr/bin/perl
#use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use StructureLibrary::Sequence;
use StructureLibrary::Tools;
use StructureLibrary::Structure;
use StructureLibrary::RNAtoolsConfig; 
use vars qw(%CONFIG);
*CONFIG = \%RNAtoolsConfig::CONFIG;


=head1 NAME
structure_accuracy.pl 

=head2 SYNOPSIS
Calculate accuracy values of a given structure for all occurrences of a given motif.

COMPULSORY
	-msfas=s		motif sequence + structure fasta file
	-str-f=s		structure format [DB,C]
	-dot=s			Dotplot output from folding algorithm
	
OPTIONAL
	-help			Present this message
	-rsfas=s		Reference sequence fasta file	
	-pu=s AND -puh=s Unpaired probability file and number of header lines in file
	-ed				Calculate the ensemble energy difference and probability of
					constraint structure in -str-f using RNAfold
	-o				Output file name
					
=head3 DESCRIPTION

You have a motif that occurs once or often within a sequence and you know the precise
structure that occurs within the motif. You have also folded this longer sequence with  
a method of choice. Then you can add all these inputs and the script will give you
different accuracy scores of the structure for each position the motif occurs in the
longer sequence. The structure and the motif must have the same length! 

NOTE: use absolute paths to your files if using the cluster to calculate your results!

AUTHOR(S): Sita Lange and Daniel Maticzka

=head4 OPTIONS

        -help   This message.
        
        -msfas	<MOTIF_STRUCTURE_FASTA> e.g. motif-structure.fasta (COMPULSORY) 
        		This file in fasta format includes the name of the motif in the
        		header, the motif sequence in the first row after the header
        		and the structure in the format given by -str-f in the 3rd row. 
        		See -str-f for more information on the format. for the allowed
        		symbols in the dot-bracket format.
        		
       	-str-f	<STRUCTURE_FORMAT> e.g. DB (COMPULSORY) 
       			This gives the format of your structure in -msfas. The options are:
       			
       			(1) DB = DOT-BRACKET = "().-"
       				matching brackets ( ): base i pairs with base j 
       				. : base i is unpaired
       				- : structure of base i is unkown
       				
       			(2) C = CONSTRAINT = "()x.|<>", which are the symbols for RNAfold
       				matching brackets ( ): base i pairs with base j
       			 	| : paired with another base
 					. : no constraint at all
 					x : base must not pair
 					< : base i is paired with a base j<i
 					> : base i is paired with a base j>i
 					
        -dot    <DOTPLOT> e.g. "seq_dp.ps" (COMPULSORY)
        		The dotplot file in ps format (e.g. output of RNAplfold).
        		
        -rsfas <REFERENCE_SEQUENCE_FASTA> e.g. mRNA.fasta (OPTIONAL) 
        		This is a fasta file that includes the sequence of your reference 
        		sequence. The reference sequence is the (longer) sequence in which
        		your structural motif occurs. 
        		NOTE: this file must inlcude only 1 sequence entry!
        		
        -pu		<UNPAIRED> e.g. "seq_pu.acc"
        		The unpaired probabilities that fit to the dotplot in -dot. This means
        		they are not specifically calculated. If RNAfold was used there are
        		no unpaired probabilities and these can be calculated by not giving
        		this option
        		
        -puh	<NUM_HEADER> e.g. 2 (COMPULSORY with -pu)
        		The number or header lines in the PU file.
        		
        -ed		Calculate the ensemble energy difference and the probability of the
        		constraint structure(s) using the constraint given in -msfas and 
        		RNAfold. NOTE: if DOTBRACKET format is given, then this will be converted
        		into the CONSTRAINT format.  
        		NOTE: The reference fasta sequence -rsfas is compulsory if you want
        		to calculate the energy differences and probabilities.
        
        -o		<FILE_NAME> e.g. struct.acc
        		Name of the output file. If this is not given, the output is written
        		to STDOUT.
        
        -sge	<Username>
        		The script is run on the cluster for faster calculation using your
        		username. 
        		
=head5 EXAMPLES



=cut


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################

# command line options
my ($i_help, $i_man, $i_ref_seq_fasta, $i_motif_str_fasta, $i_str_format, $i_dot, $i_pu, 
	$i_pu_num_header, $i_RNAfoldED, $i_RNAfold_path, $i_out, $i_sge);

my $options = GetOptions ("help"	=> \$i_help,
                         "man"  	=> \$i_man,
                         "rsfas=s"	=> \$i_ref_seq_fasta,
                         "msfas=s"	=> \$i_motif_str_fasta,
                         "str-f=s"	=> \$i_str_format,
                         "dot=s" 	=> \$i_dot,
                         "pu=s"		=> \$i_pu,
                         "puh=s"	=> \$i_pu_num_header,
                         "ed"		=> \$i_RNAfoldED,
                         "rnafold=s"=> \$i_RNAfold_path,
                         "o=s"		=> \$i_out,
                         "sge"		=> \$i_sge
							);
                      
                      
# check input
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;
pod2usage(-exitstatus => 0, -verbose => 2) if $i_man;
($options) or pod2usage(2);

# has a structure been given?
($i_motif_str_fasta) or pod2usage("INPUT ERROR: -msfas is mandatory!\n");
# has the structure format been given?
($i_str_format) or pod2usage("INPUT ERROR: -str-f is a compulsory parameter!\n");
$i_str_format = Structure::check_str_format(0, $i_str_format, "script structure_accuracy.pl");
if($i_RNAfoldED){
	($i_ref_seq_fasta) or pod2usage("INPUT ERROR: -rsfas <REFERENCE_FASTA> is compulsory\n".
					"with option -ed!\n");
}
# has the dotplot file been given?
($i_dot) or pod2usage("INPUT ERROR: -dot DOTPLOT is mandatory");
(-e $i_dot) or pod2usage("INPUT ERROR: $i_dot must be an existing file\n");

## check that if -pu is given, then -puh is also given
#if($i_pu){
#	($i_pu_num_header) or pod2usage("INPUT ERROR: -pu has been given, so -puh must be given too \n(number of header lines in unpaired probability file)!\n");
#}


		
###############################################################################
# GLOBAL VARIABLES
###############################################################################
my $CURRUSER = getlogin;
my $CURRDIR = getcwd;
my $TMPDIR = "$CURRDIR/$CONFIG{TMPDIR}";
(-e $TMPDIR) or system("mkdir $TMPDIR");
my $SGE_SUBMIT = "/home/sita/bin/RNAtools/StructureLibrary/generic_submit_to_cluster.sh";
#my $RNAFOLDPATH = "/usr/local/vrna/1.8.4/bin/RNAfold";

my $REFERENCE_SEQ_HREF;
my $MOTIFS_STRUCTURE_HREF;
my $MOTIF_ORDER;
my $NUM_MOTIFS;

#use vars qw(%STR_FORMATS $DOTBRACKET $CONSTRAINT);

###############################################################################
# EXECUTION CODE
###############################################################################

# see if this job should be run on the cluster

if($i_sge){
	system("ssh $CURRUSER\@biui.informatik.uni-freiburg.de '".
	"export SGE_ROOT=/opt/sge-6.0; ".
	"/opt/sge-6.0/bin/lx24-amd64/qsub ".
	"$SGE_SUBMIT $CURRDIR \"$options\"'");
	exit;
}
		

# contents for output file
my @out_contents = ();

# read structure file and check that the symbols of both sequence and structure
my $motif_structure_contents_aref = Structure::read_seq_structure_fasta_file(
											$i_motif_str_fasta, $i_str_format);
$MOTIFS_STRUCTURE_HREF = $motif_structure_contents_aref->[0];
$MOTIF_ORDER = $motif_structure_contents_aref->[1];
$NUM_MOTIFS = keys (%{$MOTIFS_STRUCTURE_HREF});
($NUM_MOTIFS) or die "ERROR: there are no motifs, either the parsing of the motif".
		" file failed, or the file $i_motif_str_fasta is empty!\n";

# set reference sequence and check that there is only one 
if($i_ref_seq_fasta){
	(-e $i_ref_seq_fasta) or die("INPUT ERROR: $i_ref_seq_fasta must be an existing file!\n");
	my $ref_seq_contents_aref = Sequence::read_fasta_file($i_ref_seq_fasta);
	$REFERENCE_SEQ_HREF = $ref_seq_contents_aref->[0];
	my $n = keys (%{$REFERENCE_SEQ_HREF});
	($n==1) or die
			("INPUT ERROR: There is more than one element in $i_ref_seq_fasta");
			
#	# check reference sequence
#	print STDERR $REFERENCE_SEQ_HREF->{$ref_seq_contents_aref->[1]->[0]}."\n";
}



# get base-pair probabilities
my $dot_aref = Structure::parse_dotplot_return_bp_probs_and_sequence($i_dot);
my $bp_prob_href = $dot_aref->[0];
# translate sequence in dotplot to RNA
$dot_aref->[1] =~ tr/T/U/d;

my $ol = "refseqID\tmID\tpos\tscore_type\tscore";
print STDOUT $ol."\n";
push(@out_contents, $ol);

# process each motif one-by-one
foreach my $mID (keys %{$MOTIFS_STRUCTURE_HREF}){
	
	my $motif_occ_aref;
	my $refseq 			= "NA";
	my $refseqID 		= "NA";
	
	my $mlen = length($MOTIFS_STRUCTURE_HREF->{$mID}->{SEQUENCE});
	
	#translate motif into RNA 
	$MOTIFS_STRUCTURE_HREF->{$mID}->{SEQUENCE} =~ tr/T/U/d;
	
	# find occurences of the motif in the respective sequence
	# a reference sequence is given, find the starting positions
	if($REFERENCE_SEQ_HREF){
		my @a 		= keys(%{$REFERENCE_SEQ_HREF});
		$refseqID 	= shift(@a);
		$refseq 	= $REFERENCE_SEQ_HREF->{$refseqID};
		#translate reference sequence into RNA 
		$refseq 	=~ tr/T/U/d;
		
		($refseq =~ /^[AGUC]+$/) or 
			die ("INPUT ERROR: The given sequence in $i_ref_seq_fasta does not consist".
				" of RNA (or DNA)!\n");
		
		# find occurences of motif in reference sequence
		$motif_occ_aref = Sequence::find_motif($refseq, 
			$MOTIFS_STRUCTURE_HREF->{$mID}->{SEQUENCE});
		
		($motif_occ_aref) or 
			print STDERR ("WARNING: No occurrences of $mID was found in sequence".
						 " $refseqID!\n");
		
#		# print occurences to check
#		if($motif_occ_aref){
#			my $first = shift(@{$motif_occ_aref});
#			my $regions = ($first+1)."-".($first+$mlen);
#			foreach (@$motif_occ_aref){
#				$regions .= ",".($_+1)."-".($_+$mlen);
#			}
#			print $regions."\n";
#		}
		
		die "INPUT ERROR: The reference sequence does not correspond to the sequence".
			" in the dotplot!\n".
			"REFERENCE SEQUENCE: $refseq\nDOTPLOT SEQUENCE: $dot_aref->[1]\n"
			if($dot_aref->[1] ne $refseq);
		
		
	# no reference sequence is given, the reference sequence is the motif 
	# sequence and the starting position is at 0, counting from 0 to n-1.
	# Count from 0 as this is what Sequence::find_motif() does.
	} else {
		push( @{$motif_occ_aref}, 0);
		$refseq = $MOTIFS_STRUCTURE_HREF->{$mID}->{SEQUENCE};
		
		die "INPUT ERROR: The motif sequence does not correspond to the sequence".
			" in the dotplot!\n".
			"MOTIF SEQUENCE (ALSO REFERENCE): $MOTIFS_STRUCTURE_HREF->{$mID}->{SEQUENCE}\n".
			"DOTPLOT SEQUENCE: $dot_aref->[1]\n"
			if ($dot_aref->[1] ne $refseq);
	}
	my $refseq_len = length($refseq);
	
	# each occurence is calculated separately
	my $ed_res_href;
	foreach my $occ (@{$motif_occ_aref}){
		
		$occ += 1; # now counting 1 to n and not 0 to n-1
		
		# calculate accuracy score
		my $acc = Structure::accuracy($MOTIFS_STRUCTURE_HREF->{$mID}->{STRUCTURE}, 
					$i_str_format, $occ, $bp_prob_href,0,$i_pu, $i_pu_num_header,
					$refseq, $refseq_len,"normalise",1);

					
		$ol = "$refseqID\t$mID\t$occ\tNormalised Accuracy\t$acc";
		
		print STDOUT $ol."\n";
		push(@out_contents, $ol);
		
		
		
		# calculate ensemble energies etc
		my $constraint_struct = $MOTIFS_STRUCTURE_HREF->{$mID}->{STRUCTURE};
		if($i_RNAfoldED){
			if($i_str_format eq $CONFIG{STR_FORMATS}->{DOTBRACKET}){
				print STDERR "WARNING: converting structure format from dotbracket to constraint (RNAfold -C)!\n";
				print "before conversion: $MOTIFS_STRUCTURE_HREF->{$mID}->{STRUCTURE}\n";
				$constraint_struct = Structure::convert_structure_format($MOTIFS_STRUCTURE_HREF->{$mID}->{STRUCTURE},
					$CONFIG{STR_FORMATS}->{DOTBRACKET}, $CONFIG{STR_FORMATS}->{CONSTRAINT});
				
				$i_str_format = $CONFIG{STR_FORMATS}->{CONSTRAINT};
				print "after conversion: ".$constraint_struct."\n";
			}
			
			if($i_str_format eq $CONFIG{STR_FORMATS}->{CONSTRAINT}){
				# create constraint fasta file for this motif occurence
				my $tempfile = "$TMPDIR/temp.constraint.fasta";
				my @cons_fasta_a = (); 
				push(@cons_fasta_a, "> ${mID}_".${occ});
				push(@cons_fasta_a, $refseq);
				
				$constraint_struct = create_x_dots($occ-1).
				$constraint_struct.
				create_x_dots($refseq_len - ($occ + $mlen - 1));
				push(@cons_fasta_a, $constraint_struct);
				
				Tools::write_file_from_array($tempfile, \@cons_fasta_a);
				
				# call method to calculate RNAfold results
				$ed_res_href = ensemble_energies_and_probability_with_constraint($i_ref_seq_fasta,
					$tempfile, $CONFIG{RNAFOLDPARAMS}, $TMPDIR,
					1, 1, $ed_res_href);
					
				$ol = "$refseqID\t$mID\t$occ\tED\t$ed_res_href->{ED}";
				print STDOUT $ol."\n";
				push(@out_contents, $ol);
				$ol = "$refseqID\t$mID\t$occ\tConstraint Probability\t".
						"$ed_res_href->{CONSTRAINT_P}";
				print STDOUT $ol."\n";
				push(@out_contents, $ol);
			}
		}
	}
	
	# write output file, if wanted
	if($i_out){
		Tools::write_file_from_array($i_out, \@out_contents);
	}
	
	# remove all temporary files
	system ("rm -R $TMPDIR/*");
}



###############################################################################
# METHODS
###############################################################################

# create a string with x dots, e.g. '....' has four dots
sub create_x_dots{
	my $x = shift;
	
	my $string = "";
	
	for(my $i=0 ; $i < $x ; $i++){
		$string .= ".";
	}
	return $string;
}