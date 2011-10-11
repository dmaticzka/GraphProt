#!/usr/bin/perl
#use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;


=head1 NAME

Dotplot Viewer

=head1 SYNOPSIS

dotview.pl -dot DOTPLOT [OPTIONS]


    This script takes the output dotplot in ps format from the Vienna package
    programs and highlights interesting regions in it. It is possible to
    mark subsequences, regions given by sequence positions, and a distinct
    structure of interest. The subsequences of interest are written in capital
    letters, the rest in lower case. You can also add a structure independent of the
    structure in the dotplot as outlined dots to see how it compares with the
    dotplot. Furthermore, it is possible to highlight the dots in different
    colours, according to what you define, default is red. Finally, as
    dotplots tend to get too big for viewers for longer sequences, it is
    possible to crop the dotplot to a region you are interested in. This
    retains the same structure, but shows only a subplot, similar to zooming
    in, but in a static way for one defined region.

NOTE: The highlighting only works if you enter the dotplot in the original format as 
produced by RNAplfold.

AUTHOR: Sita Lange

OPTIONS

		
        -help   This message.
        
        -dot    <DOTPLOT> e.g. "seq_dp.ps"
        		The dotplot file in ps format (output of RNAplfold).
        
        -o		<FILE_NAME>
        		Name of the output file.
        
        -seq	[SUBSEQUENCES] e.g. "ATTC,GGGC,GUGUC"
        		A list of subsequences separated by commas.
        
        -reg	[STARTPOS-STOPPOS] e.g. "1-3,5-7,9-11"
        		A list of regions where all dots involved in region are 
        		highlighted according to -col, given as start and stop positions, 
        
        		If -reg and -seq are given at the same time, then regions and 
        		subsequences are considered separately and the subsequences
        		are only highlighted by capital letters and the dots within the regions
        		are marked with the given colours, without marking the corresponding
        		subsequences. This means subsequences of interest can be highlighted 
        		independently of local structures of interest. 		
        
        -col	[COLOURS] e.g. "red,green,blue"
        		A list of colours according to which the subsequence or the region lists 
        		shall be highlighted. If this list is shorter to the the aforementioned 
        		list, the colours will be repeated in the given order.
        		Possible colours include: red, green, blue, lightblue, 
        		magenta, yellow, orange, purple, darkgreen.
        		Default is to highlight every region red.
        
        -crop	<STARTPOS>-<STOPPOS> e.g. "20-30"
        		Crop the dotplot to given start and stop position, e.g. -crop 50-100 crops 
        		the given dotplot to only structures within the subsequence between 50  
        		to 100, including these positions. 
        
        -str	(<COLOUR>, <TYPE>, [BASEPAIRS]) e.g. "green,1:10,2:9,3:8,4:7"
        		The colour, followed by a list of base pairs that make up an exact 
        		structure of interest, which shall be highlighted by the given colour. 
        		Base 1 pairing with base 10 is given as such, 1:10 and only base-pairs that 
        		exist are highlighted, the others are ignored.
        		
        -str-db (<COLOUR>,<START>,<DOT BRACKET STRUCTURE>) 
        			e.g. "purple, 20, (((((...)))..((...)).))"
        		A structure given in dot-bracket notation (only round brackets and dots).
        		This can be used instead of giving the exact base pairs in -str. You
        		also need to give the colour and the starting position of the structure.
        		
        
        -ne-str Plot non-existant structure: It is possible to set this flag if you want to  
        		add a non-existant structure into the dotplot. The structure given in -str 
        		is then drawn into the dotplot in the given colour but as unfilled 
        		rectangles of a fixed  with an existing structure, the original dots remain 
        		visible in the background. size. Only valid in combination with -str
        		and -str-db.
        		
        		
        -strict If this flag is given, only base pairs with both bases within the given
        		regions are coloured. The default would be to highlight all base-pairs
        		with at least one base falling within the region.
        		
        -pdf 	The output of the resulting dotplot is a PDF
        
        -png	The output of the resulting dotplot is a PNG
        
        -r		The original ps file is replaced by the new dotplot. This means the 
        		original file given in -dot is deleted.
        		
        		
NOTE: all indices on the sequence are given from 1..n and not from 0..(n-1), 
		where n is the length of the sequence.
        	
        		
EXAMPLES
	
	(1) dotview.pl -d seq_dp.ps -seq "ATTC,GGGC,GUGUC" 
		
		Highlights the dots within all occurrences of the given subsequences red 
		(default) and writes the subsequences of interest in capital letters, 
		whereas the rest are in lower case.
	
	
	(2) dotview.pl -d seq_dp.ps -seq "ATTC,GGGC,GUGUC" -crop "20-50" -pdf -r
		
		Does the same as in (1) and crops the dot plot to positions 20 to 50 in the 
		original sequences. This is useful if the original dotplot is too large, 
		i.e. it is a type of zoom function. The resulting dotplot is produced in 
		PDF format. Finally the original file seq_dp.ps is deleted.
	
	
	(3) dotview.pl -d seq_dp.ps -reg "1-3,5-7,9-11" -col "red,green,blue"
		
		Colours all dots in region 1-3 red, 5-7 green, and 9-11 blue.
		If the regions overlap, then a dot is given the colour of the last region.
		Also marks the subsequences within the given regions with capital letters.
		
		
	(4) dotview.pl -d seq_dp.ps -reg "1-3,5-7,9-11" -col "red,green,blue" 
			-seq "ATTC,GGGC,GUGUC" 
		
	  	Does the same as (3) but does not mark the regions with capital letters, but
	  	instead marks the given subsequences with capital letters. The dots within
	  	the subsequences are then not highlighed by a different colour.
	  	
	 (5) dotview.pl -d seq_dp.ps -str-db "blue,10,(((....)))" -ne-str
	 
	 	Plots an additional, non-existant structure into the dotplot, with blue dots
	 	that are only outlines and not filled like the real structure. The structure
	 	corresponds to the base-pairs 10:19,11:18,12:17.

=head1 DESCRIPTION

=cut


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################

# command line options
my ($i_help, $i_man, $i_dot, $i_out, $i_seq, $i_reg, $i_col, $i_crop, $i_str, 
		$i_str_db, $i_ne_str, $i_strict, $i_pdf, $i_png, $i_replace);

my $options = GetOptions ("help"	=> \$i_help,
                         "man"  	=> \$i_man,
                         "dot=s" 	=> \$i_dot,
                         "o=s"		=> \$i_out,
                         "seq=s"	=> \$i_seq,
                         "reg=s"	=> \$i_reg,
                         "col=s"	=> \$i_col,
                         "crop=s"	=> \$i_crop,
                         "str=s"	=> \$i_str,
                         "str-db=s" => \$i_str_db,
                         "strict"	=> \$i_strict,
                         "ne-str"	=> \$i_ne_str,
                         "pdf"		=> \$i_pdf,
                         "png"		=> \$i_png,
                         "r"		=> \$i_replace);
                      
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;
pod2usage(-exitstatus => 0, -verbose => 2) if $i_man;
($options) or pod2usage(2);
($i_dot) or pod2usage("INPUT ERROR: -dot DOTPLOT is mandatory");
(-e $i_dot) or pod2usage("INPUT ERROR: $i_dot must be an existing file\n");
($i_seq ||$i_reg || $i_crop || $i_str || $i_str_db) or pod2usage("INPUT ERROR: at least one of the following options".
		" should be given: -seq [SUBSEQUENCES], -reg [REGIONS], -crop STARTPOS-STOPPOS, ".
		"-str (COLOUR,[BASEPAIRS]), -str-db (<COLOUR>,<START>,<DOT BRACKET STRUCTURE>)");
pod2usage("INPUT ERROR: only one of the following options can be given: -str or -str-db!\n") 
			if ($i_str && $i_str_db);
			
pod2usage("INPUT ERROR: only one of the following options can be given: -pdf or -png!\n") 
			if ($i_pdf && $i_png);
		
# delete spaces
$i_reg =~ tr/ //d if($i_reg);
$i_seq =~ tr/ //d if($i_seq);
$i_crop =~ tr/ //d if($i_crop);
$i_col =~ tr/ //d if($i_col);
$i_str =~ tr/ //d if($i_str);
$i_str_db =~ tr/ //d if ($i_str_db);

# translate given sequence to RNA
$i_seq = uc($i_seq);
$i_seq =~ tr/T/U/d if ($i_seq);

		
###############################################################################
# GLOBAL VARIABLES
###############################################################################
my $TMPDIR = "/scratch/1/sita/tmp/";

#ps colour definitions
my $RED 		= "1 0 0 setrgbcolor";	
my $GREEN		= "0 1 0 setrgbcolor";
my $BLUE	= "0 0 1 setrgbcolor";
my $LIGHTBLUE	= "0 1 1 setrgbcolor";
my $MAGENTA		= "1 0 1 setrgbcolor";
my $YELLOW		= "1 1 0 setrgbcolor";
my $WHITE		= "1 1 1 setrgbcolor";
my $BLACK		= "0 0 0 setrgbcolor";
my $ORANGE		= "1.0 0.7 0.0 setrgbcolor";
my $PURPLE		= "0.7 0.3 1.0 setrgbcolor";
my $BROWN		= "0.7 0.3 0.0 setrgbcolor";
my $DARKGREEN	= "0.0 0.5 0.0 setrgbcolor";
my $DEFAULTCOLOUR = $RED;

my $SEQUENCE = "";
my $SEQ_BOUND_LEFT = 1;
my $SEQ_BOUND_RIGHT = 1;

my $PS_RECT_MACRO = 	"/rect { %size x y rect - draws rectangle centered on x,y\n".
  		 				"2 index 0.5 mul sub            % x -= 0.5\n".
   						"exch 2 index 0.5 mul sub exch  % y -= 0.5\n".
   						"3 -1 roll dup rectstroke\n".
						"} bind def";


my $PS_OBOX_MACRO = 	"/obox {\n".
   						"logscale {\n".
      					"log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n".
   						"} if\n".
   						"3 1 roll\n".
   						"exch len exch sub 1 add rect\n".
						"} bind def";
						
my $PS_LOBOX_MACRO= "/lobox {".
   					 "3 1 roll\n".
   					 "len exch sub 1 add rect\n".
					 "} bind def";


my $PS_BOX_OUTLINED = "obox";
my $PS_BOX_FILLED = "ubox";
my $PS_BOX_MFE = "lbox";
my $PS_BOX_MFE_OUTLINED = "lobox";


my $GLOBALPS = 0;

###############################################################################
# EXECUTION CODE
###############################################################################

dotbracket2basepairs() if ($i_str_db);

my $ps_file_contents_aref = convert_dotplot();

# write file and do cleaning up as given by input parameters
if($i_dot =~ /(.+)\.ps/){
	my $prefix = $1.".convert";
	my $psoutput = $prefix.".ps";
	if($i_out){
		if($i_out =~ /(.+)\.p[sngdf]+$/){
			$prefix = $1;
			$psoutput = $prefix.".ps";
		} else {
			$prefix = $i_out;
			$psoutput = $i_out.".ps";
		}
	}
	# write ps file	
	write_file_from_array($psoutput, $ps_file_contents_aref);
	# replace ps with pdf
	if($i_pdf){
		system("ps2pdf -dEPSCrop $psoutput $prefix.pdf");
		system("rm $psoutput");
	
	# replace ps with png
	} elsif($i_png) {
		if($i_out && ($i_out =~ /\.png$/)){
			system("convert -density 1000 $psoutput $i_out");
		} else {
			system("convert -density 1000 $psoutput $prefix.png");
		}
	}
	# delete original dotplot
	if($i_replace){
		system("rm $i_dot");
	}
} else {
	pod2usage("INPUT ERROR: The dotplot file $i_dot is not a ps file, i.e. it does not end in '.ps'!\n");
}
print STDERR "Finished!\n";

###############################################################################
# METHODS
###############################################################################

###############################################################################
# Here the main convertion of the dotplot takes place. Every input parameter
# is controlled and carried out.
# INPUT: none
# OUTPUT: Array of contents for the new PS file.
###############################################################################
sub convert_dotplot{
	
	# converted file contents in array	
	my (@filecontents) = ();
	my $file_ending = "";
	
	# create array of colours in ps format
	my @col_array = ();
	if($i_reg || $i_seq){ # only if regions or subsequences are given
		if($i_col){
			@col_array = split(',', $i_col);
			my $n = @col_array;
			for(my $i = 0 ; $i < $n ; $i++){
				$col_array[$i] = convert_colours_to_ps_string($col_array[$i]);
			}
		} else {
			@col_array = ($DEFAULTCOLOUR); # set default
		}
	}
	
	# base-pairs B(i,j) to be highlighted with given colours saved in a hash of hashes
	# in the form of $h{i}{j} = colour, with i < j.
	my $struct_colours_hhref = {};
	my @i_reg_array = ();
	if($i_reg){
		@i_reg_array = split(',', $i_reg);
		create_structure_colour_hash($struct_colours_hhref, \@i_reg_array, \@col_array);
	}
	# get colours for base pairs in structure, if $i_ne_str is not given.
	if($i_str && !$i_ne_str){
		my @i_struct_array = split(",", $i_str);
		my $col = shift(@i_struct_array);
		pod2usage("INPUT ERROR: The structure in -str has been given in the wrong format,
			there is no leading colour: col=$col in $i_str!\n") if($col =~ /\d+/);
		my @cola = (convert_colours_to_ps_string($col));
		create_structure_colour_hash($struct_colours_hhref, \@i_struct_array, \@cola);
	} 
	
	# create subsequences array
	my @i_seq_array = ();
	if ($i_seq){
		@i_seq_array = split(",", $i_seq);
	}
	
	# open ps file for parsing
	open(IN_HANDLE, "<$i_dot") || die "couldn't open the following file: $i_dot/n";
									 
		my ($line, $seqline) = "";	
		my $inseq = 0;
		while ($line = <IN_HANDLE>){
			chomp ($line);
			
			# beginning of macro definitions for boxes, add extra definitions
			if($line =~ /\/lpmin/){
				push(@filecontents, $PS_RECT_MACRO);
				push(@filecontents, "");
				push(@filecontents, $PS_OBOX_MACRO);
				push(@filecontents, "");
				push(@filecontents, $PS_LOBOX_MACRO);
			}
			
			# get sequence and convert it according to the given parameters
			if($line =~ /\/sequence/){
				$inseq = 1;
				$seqline = $line;
			} elsif($inseq){ # in sequence region
				if($line =~ /def/){ # end of sequence region
					$inseq = 0;
					$seqline .= $line;
					
					set_sequence($seqline);
					
					# if subsequences are given, but no regions are given, highlight base-pairs 
					# within subsequence regions
					if($i_seq){
						my $seqregions_aaref = convert_subsequences_to_regions(\@i_seq_array);
						push(@filecontents,mark_and_crop_sequence($seqregions_aaref));
						create_structure_colour_hash($struct_colours_hhref, $seqregions_aaref, 
						\@col_array) if (!$i_reg);
					
					# no subsequences are given, but subsequences indicated by the regions need to be marked.
					} elsif($i_reg) {
						push(@filecontents,mark_and_crop_sequence(\@i_reg_array));
					
					# no regions to be marked, but the sequence needs to be cropped
					} else {
						# method marks the regions in upper case and also sets the 
						# originial sequence and the boundaries given by crop in the corresponding
						# global variables.
						push(@filecontents,mark_and_crop_sequence());  
					}
				} else {
					$seqline .= $line;
				}
			
			# get the base pairs, crop and colour them according to the given parameters
			} elsif($line =~ /\d+\s\d+\s\d(\.?\d*)\s$PS_BOX_FILLED/){
				
				my @line_a = split (" ", $line); 
				# the originial file should have 4 columns: base i, base j, sqrt of bp prob, ubox
				my $c = @line_a;
				pod2usage( "INPUT ERROR: The input file $i_dot is not an original dotplot file\n")
						if($c != 4);
				my $basei = $line_a[0];
				my $basej = $line_a[1];
				my $prob  = $line_a[2];
				if($basei > $basej){ #only if the second base is larger than the first
					my $temp = $basei;
					$basei = $basej;
					$basej = $temp;
				}
				
				# colour and crop base-pairs
				if($basei >= $SEQ_BOUND_LEFT && $basej <= $SEQ_BOUND_RIGHT){
					
					# modify basei and base j to cropped indexes, if applicable
					my $ori_basei = $basei;
					my $ori_basej = $basej;
					if($i_crop){
						$basei = $ori_basei - $SEQ_BOUND_LEFT + 1;
						$basej = $ori_basej - $SEQ_BOUND_LEFT + 1;
					}
					
					# highlighted base-pair
					if(defined($struct_colours_hhref->{$ori_basei}->{$ori_basej})){
						push(@filecontents, $struct_colours_hhref->{$ori_basei}->{$ori_basej}.
						" $basei $basej $prob $PS_BOX_FILLED");
					
					# highlighted base i
					} elsif (defined($struct_colours_hhref->{$ori_basei}->{"N"})){
						push(@filecontents, $struct_colours_hhref->{$ori_basei}->{"N"}.
						" $basei $basej $prob $PS_BOX_FILLED");
						
					# highlighted base j
					} elsif (defined($struct_colours_hhref->{$ori_basej}->{"N"})){
						push(@filecontents, $struct_colours_hhref->{$ori_basej}->{"N"}.
						" $basei $basej $prob $PS_BOX_FILLED");
						
					# non-highlighted base-pair
					} else {
						push(@filecontents, $BLACK." $basei $basej $prob $PS_BOX_FILLED");
					}
				}
				
			# IF global dotplot, then there lbox is defined for the MFE structure
			# In this case these are given the same markings as the probabilities,
			# but a different handling is possible here.
			} elsif($line =~ /\d+\s\d+\s\d\.\d+\s$PS_BOX_MFE/){
				
				# set global ps setting
				$GLOBALPS = 1;
				
				my @line_a = split (" ", $line); 
				# the originial file should have 4 columns: base i, base j, sqrt of bp prob, ubox
				my $c = @line_a;
				pod2usage( "INPUT ERROR: The input file $i_dot is not an original dotplot file\n")
						if($c != 4);
				my $basei = $line_a[0];
				my $basej = $line_a[1];
				my $prob  = $line_a[2];
				if($basei > $basej){ #only if the second base is larger than the first
					my $temp = $basei;
					$basei = $basej;
					$basej = $temp;
				}
				
				# colour and crop base-pairs
				if($basei >= $SEQ_BOUND_LEFT && $basej <= $SEQ_BOUND_RIGHT){
					
					# modify basei and base j to cropped indexes, if applicable
					my $ori_basei = $basei;
					my $ori_basej = $basej;
					if($i_crop){
						$basei = $ori_basei - $SEQ_BOUND_LEFT + 1;
						$basej = $ori_basej - $SEQ_BOUND_LEFT + 1;
					}
					
					# highlighted base-pair
					if(defined($struct_colours_hhref->{$ori_basei}->{$ori_basej})){
						push(@filecontents, $struct_colours_hhref->{$ori_basei}->{$ori_basej}.
						" $basei $basej $prob $PS_BOX_MFE");
					
					# highlighted base
					} elsif (defined($struct_colours_hhref->{$ori_basei}->{"N"})){
						push(@filecontents, $struct_colours_hhref->{$ori_basei}->{"N"}.
						" $basei $basej $prob $PS_BOX_MFE");
						
					} elsif (defined($struct_colours_hhref->{$ori_basej}->{"N"})){
						push(@filecontents, $struct_colours_hhref->{$ori_basej}->{"N"}.
						" $basei $basej $prob $PS_BOX_MFE");
						
					# non-highlighted base-pair
					} else {
						push(@filecontents, $BLACK." $basei $basej $prob $PS_BOX_MFE");
					}
				}
				
				
			} elsif($line =~ /showpage/){
				$file_ending = $line;
			} elsif($line =~ /end/){
				$file_ending .= "\n".$line;
			} elsif($line =~ /EOF/){
				$file_ending .= "\n".$line;
			} else {
				push(@filecontents, $line);
			}
		}
		
	# add non-existing structure into the plot
	# get colours for base pairs in structure, if $i_ne_str is  given.
	if($i_str && $i_ne_str){
		my @i_struct_array = split(",", $i_str);
		my $col = shift(@i_struct_array);
		my @cola = (convert_colours_to_ps_string($col));
		my $ne_struct_colours_hhref = "";
		create_structure_colour_hash($ne_struct_colours_hhref, \@i_struct_array, \@cola);
		
		my @keys = keys (%{$ne_struct_colours_hhref});
		my $num = @keys;
		
		foreach my $bp (@i_struct_array){
			if($bp =~ /(\d+):(\d+)/){
				my $basei = $1;
				my $basej = $2;
				if($basei > $basej){
					my $tmp = $basei;
					$basei = $basej;
					$basej = $tmp;
				}
				if($i_crop){
					if($basei >= $SEQ_BOUND_LEFT && $basej <= $SEQ_BOUND_RIGHT){
						push(@filecontents, $ne_struct_colours_hhref->{$basei}->{$basej}.
						" ".($basei-$SEQ_BOUND_LEFT+1)." ".($basej-$SEQ_BOUND_LEFT+1)." 1 $PS_BOX_OUTLINED");
						push(@filecontents, $ne_struct_colours_hhref->{$basei}->{$basej}.
						" ".($basei-$SEQ_BOUND_LEFT+1)." ".($basej-$SEQ_BOUND_LEFT+1)." 1 $PS_BOX_MFE_OUTLINED") if ($GLOBALPS);
					}
				} else {
					push(@filecontents, $ne_struct_colours_hhref->{$basei}->{$basej}.
					" ".$basei." ".$basej." 1 $PS_BOX_OUTLINED");
					push(@filecontents, $ne_struct_colours_hhref->{$basei}->{$basej}.
					" ".$basei." ".$basej." 1 $PS_BOX_MFE_OUTLINED") if ($GLOBALPS);
				}
			} else {
				pod2usage("Wrong format for base pairs: base-pair $bp");
			}
		}

	} 
	# push end of file into contents and wind up method
	push(@filecontents, $file_ending);
	close(IN_HANDLE);
	return \@filecontents;
}

#######################################################################
# Gets a list of regions and marks these within the sequence as 
# upper case letters. The remaining letters are written in lower
# case. It also crops the sequence to the given index range in -crop
# and transforms the sequence back into ps format with a line
# break every 60 nucleotides.
# INPUT: regions_aref The array of regions
# OUTPUT: the marked and cropped sequence in PS format
#######################################################################
sub mark_and_crop_sequence{
	my ($regions_aref) = @_;
	
	my $marked_seq = lc($SEQUENCE);
	
	# set borders
	if($i_crop){
		my @crop_array = split("-", $i_crop);
		pod2usage("INPUT ERROR: Can not perform -crop $i_crop. Region is either \n".
				"out of bounds or backwards.\n")
			if(($crop_array[0] >= $crop_array[1]) || ($crop_array[0] < 1) || 
				($crop_array[1] > $SEQ_BOUND_RIGHT));
		$SEQ_BOUND_LEFT = $crop_array[0];
		$SEQ_BOUND_RIGHT = $crop_array[1];
	}
	
	# mark the sequence
	foreach my $region (@{$regions_aref}){
		
		my @reg_array = split("-", $region);
		my $lb = $reg_array[0]-1;
		my $len = $reg_array[1]-$reg_array[0]+1;
		my $ss = substr($marked_seq, $lb, $len);
		my $ucss = uc($ss);
		$marked_seq =~ s/$ss/$ucss/g;
	}
	
	# crop the sequence
	if($i_crop){
		my $len = $SEQ_BOUND_RIGHT - $SEQ_BOUND_LEFT +1;
		$marked_seq = substr($marked_seq, ($SEQ_BOUND_LEFT - 1), $len);
	}
	
	#print cropped and marked sequence for reference
#	print "Subsequence $SEQ_BOUND_LEFT - $SEQ_BOUND_RIGHT: $marked_seq\n";
	
	# return the sequence in ps format
	my $ps_seq = "/sequence { (\\\n";
	
	my $ms_len = length($marked_seq);
	my $i = 0;
	
	# cut the sequence into 60nt chunks
	while ($i < $ms_len){
		my $nt_perline = 60;
		my $tmp_seq = "";
		if($i+60-1 < $ms_len){
			$tmp_seq = substr($marked_seq, $i, $nt_perline);

		} else {
			$tmp_seq = substr($marked_seq, $i);
		}
		$ps_seq .= $tmp_seq."\\\n";
		$i += $nt_perline;
	}
	$ps_seq .= ") } def";
	return $ps_seq;
}

#################################################################
# Gets the sequence in PS format and removes all the formatting
# then sets the sequence in the global variable $SEQUENCE and
# the length of the sequence in $SEQ_BOUND_RIGHT.
# INPUT: seqline The sequence line in ps format
# OUTPUT: none
##################################################################
sub set_sequence{
	my $seqline = shift;
	
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
	$seqline =~ tr/T/U/d;
	$seqline =~ tr/t/u/d;
	# set sequence information
	$SEQUENCE = $seqline;
	$SEQ_BOUND_RIGHT = length($SEQUENCE);
}

############################################################
# Just converts the given subsequences into regions, i.e. 
# index ranges on the original sequence, counting from 1..n
# INPUT subseq_aref The array of subsequences
# OUTPUT The array of regions
############################################################
sub convert_subsequences_to_regions{
	my ($subseq_aref) = @_;
	
	my $curr_index = 0;
	my $index = 0;
	my @regions = ();
	
	# iterate over all subsequences
	foreach my $subseq (@{$subseq_aref}){
		my $subseq_len = length($subseq);
		# find all occurrences of the subsequences
		while($curr_index <= $SEQ_BOUND_RIGHT){
		
			$index = index($SEQUENCE, uc($subseq), $curr_index);
			if($index == -1){
				$curr_index = $SEQ_BOUND_RIGHT+1;
			} else {
				push(@regions, ($index+1)."-".($index+$subseq_len));
				$curr_index = $index +2; # set current index to the right of occurrence
			}
		}
	}
	return \@regions;	
}

###########################################################################
# This is a simple method that takes an input colour in text form and
# converts it to the ps code for that colour.
# INPUT: colour 	The colour as an English word.
# OUTPUT: The colour in ps code.
###########################################################################
sub convert_colours_to_ps_string{
	my $colour = shift;
	
	return $RED if(uc($colour) eq "RED");
	return $GREEN if(uc($colour) eq "GREEN");
	return $BLUE if (uc($colour) eq "BLUE");
	return $LIGHTBLUE if (uc($colour) eq "LIGHTBLUE");
	return $MAGENTA if (uc($colour) eq "MAGENTA");
	return $YELLOW if (uc($colour) eq "YELLOW");
	return $BLACK if (uc($colour) eq "BLACK");
	return $ORANGE if (uc($colour) eq "ORANGE");
	return $PURPLE if (uc($colour) eq "PURPLE");
	return $BROWN if (uc($colour) eq "BROWN");
	return $DARKGREEN if (uc($colour) eq "DARKGREEN");
	
	print STDERR "WARNING: $colour is undefined, using the default colour, $DEFAULTCOLOUR\.n";
	return $DEFAULTCOLOUR;		
}

################################################################################
# The method converts the region or structure information in $str_or_reg_aref
# into colours for exact base-pairs. It differentiates between a region and
# a structure and if a region is given, then these are converted to base-pairs
# within the region. If strict is given, only base-pairs with both bases
# within the regions are marked and if not, then a base is marked and the second
# base is given as "N", which means it can indicate any other base.
# The colours for each marked base-pair are saved in a hash in the form of 
# $h{$i}{$j} = "colour code". If the colours array given is not as long
# as the regions array, then the regions/structures are coloured in permutations
# of $colours_aref.
# INPUT: 
# struct_colours_hhref 	The hash reference for the bp colours
# str_or_reg_aref 		The array ref of regions or base-pairs
# colours_aref			The array ref of colours for the regions or structures.
# OUTPUT: none
#################################################################################
sub create_structure_colour_hash{
	
	my ($struct_col_hhref, $str_or_reg_aref, $colours_aref) = @_;
	
	my $ncolours = @{$colours_aref};
	my $curr_col_i = 0; 

	# iterate through each region or base-pair
	foreach my $element (@{$str_or_reg_aref}){
		
		# the element is a region
		if($element =~ /(\d+)-(\d+)/){
			# when defining a region k-l, k<=l.
			if($1 <= $2){
				# if strict is given, only highlight base-pairs within the region
				if($i_strict){
					for(my $i = $1 ; $i < $2 ; $i++){
						for(my $j = $i+3; $j <= $2 ; $j++){
							print STDERR "WARNING: The base-pair $i:$j has been marked more than once\n" 
										if (defined ($struct_col_hhref->{$i}->{$j})&&(!$i_ne_str));
							$struct_col_hhref->{$i}->{$j} = $colours_aref->[$curr_col_i];
						}
					}
				# when strict isn't given, highlight any base pairs that have one base within the region
				# this means a base i is marked and the second base is given as "N", which means 
				# i can pair with any other base
				} else {
					for (my $i = $1; $i <=$2 ; $i++){
						print STDERR "WARNING: The base $i has been marked more than once\n"
										if (defined ($struct_col_hhref->{$i}->{"N"}));
						$struct_col_hhref->{$i}->{"N"} = $colours_aref->[$curr_col_i];
					}
				}
			} else {
				pod2usage("INPUT ERROR: In the region given by i-j, i <= j\n");
			}
		
		# the element is a base pair	
		}elsif($element =~ /(\d+):(\d+)/){
			if($1 < $2){
				print STDERR "WARNING: The base-pair $1:$2 has been marked more than once\n" 
										if (defined ($struct_col_hhref->{$1}->{$2}));
				$struct_col_hhref->{$1}->{$2} = $colours_aref->[$curr_col_i];
			} else {
				print STDERR "WARNING: The base-pair $2:$1 has been marked more than once\n" 
										if (defined ($struct_col_hhref->{$2}->{$1}));
				$struct_col_hhref->{$2}->{$1} = $colours_aref->[$curr_col_i];
			}
			
		# the element has the wrong format
		} else {
			pod2usage( "WRONG INPUT: This element from the input has an unkown format: $element\n");
		}
		
		#update colour index
		++$curr_col_i;
		$curr_col_i = 0 if($curr_col_i >= $ncolours);
			
	}
}

################################################################################
# This method parses the input given in -str-db and converts it into
# the base pair format for -str and sets the result to this variable. After
# this, the input is treated the same as if the user gave the -str input.
# INPUT: none
# OUTPUT: none
################################################################################
sub dotbracket2basepairs{
	
	my $str_basepair = "";
	
	# check format of dot-bracket structure
	if($i_str_db =~ /^(\D+),(\d+),([\(\)\.]+)$/){
		my $col = $1;
		my $start = $2;
		my @db = split("", $3);
		my $n = @db;
		my (@open, @close) = ();
		$str_basepair = $col;
		
		# parse dot-bracket notation
		for (my $i = 0; $i < $n; $i++){
			if($db[$i] eq "("){
				push(@open, ($i+$start));
			} elsif($db[$i] eq ")"){
				if(@open){
					$str_basepair .= ",".(pop(@open)).":".($i+$start) ;
				} else {
					pod2usage("INPUT ERROR: The option -str-bp has the wrong format, i.e. \n".
					"there are an uneven number of opening and closing brackets: $i_str_db!\n");
				}
			} elsif ($db[$i] ne "."){
				pod2usage("INPUT ERROR: The option -str-bp has the wrong format,\n".
						"i.e. unknown symbol $db[$i]: $i_str_db!\n");
			}
		}
		pod2usage("INPUT ERROR: The option -str-bp has the wrong format, i.e. \n".
					"there are an uneven number of opening and closing brackets: $i_str_db!\n") if(@open);
	} else {
		pod2usage("INPUT ERROR: The option -str-bp has the wrong format: $i_str_db!\n");
	}
	# set structure string 
	$i_str = $str_basepair;
}
	
##########################################################################################
# This method takes an array, which contains a string for each element that represents
# one line in the output file and writes the information in the array to  a file with the 
# given name.
# Input: file name and an array of strings
##########################################################################################
sub write_file_from_array{
	my($filename, $output_array) = @_;
	
	open(OUT_HANDLE, ">$filename") || die "couldn't open the file for writing in package Tools,".
								"sub write_file_from_array: $filename\n";
		
	chomp(@{$output_array});						
	foreach my $line (@{$output_array}){
		
		print OUT_HANDLE $line."\n";
	}
	
	close(OUT_HANDLE);
	return 1;
}
