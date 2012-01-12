#!/usr/bin/perl
use File::Basename;
use lib dirname($0);	# search skript directory for modules
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;
use StructureLibrary::Sequence;
use StructureLibrary::Structure;

=head1 NAME

createGraph.pl

=head1 SYNOPSIS

createGraph.pl -fasta=FASTA

takes sequences of fasta file and generates graph output for
fabrizio. context accessibilities are computed using a modified version
of RNAplfold that prints out a table of the context accessibilities for u=1

Options:

    -nocontext  only use accessibilities, no context probabilities
	-cutoff     do not use probabilities below or at this value (default: 0)
    -nostruct   do not compute structure part
	-fasta		fasta file including input sequences
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# get filename parts of input files
###############################################################################
#my ( $prefix, $path, $suffix ) = fileparse( $infname, "\.[^.]*" );

###############################################################################
# create temporary directory
###############################################################################
#my $tmp_template = 'template-XXXXXX';
#my $tmp_prefix = '/var/tmp/';
#my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);

###############################################################################
# readAccessibilities [DEPRECATED]
# in: id of acc file, path to files, suffix
# out: reference to array of accessibilities
###############################################################################
sub readAccessibilities {
    my ($id, $path, $suffix, $u) = @_;
    # read accessibilities
    my $fname = "$path${id}_RNAlfold.acc$suffix";
	(-f $fname) or die ("error, file does not exist: '$fname');");
	return parse_PUfile_get1u($fname, $u, 2);
}

###############################################################################
# computeAccessibilities
# in: RNA sequence
# out: reference to array of accessibilities [P,E,H,I,M]
###############################################################################
sub computeAccessibilities {
	my ($seq) = @_;
	chomp $seq;
	my $len = length($seq);
	
	# compute accessibilities
	open ACCS, "echo $seq | /home/maticzkd/src/ViennaRNA-1.8.4-context_stdout/Progs/RNAplfold " .
		"-u 1 -W $len |";
	my @accs = <ACCS>;
	close ACCS;
	
	# parse accs
	my @P;
	my @E;
	my @H;
	my @I;
	my @M;
	my $parsedseq = $accs[0];
	chomp $parsedseq;
	if ($parsedseq ne $seq) {
		die "error parsing accessibilities: expecting sequence '$seq', got '$parsedseq'"
	};
	for (my $i=1; $i<=$len; $i++) {
		my ($pos, $P, $E, $H, $I, $M) = split("\t", $accs[$i]);
		if ($pos != $i) {
			die "error parsing accessibilities: expecting position '$i', got '".$accs[$i]."'";
		}
		push @P, $P;
		push @E, $E;
		push @H, $H;
		push @I, $I;
		push @M, $M;
	}
	
	return [\@P, \@E, \@H, \@I,  \@M];
}

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $fasta;
my $path;
my $nostruct;
my $cutoff;
my $nocontext;
my $result = GetOptions (	"help"	=> \$help,
							"man"	=> \$man,
							"fasta=s" => \$fasta,
							"nostruct" => \$nostruct,
							"cutoff=f" => \$cutoff,
							"nocontext" => \$nocontext);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);
(defined $cutoff) or $cutoff = 0;

###############################################################################
# main
###############################################################################
($fasta) or die "error: specify fasta file";
(-f $fasta) or die "error: no such file '$fasta'";

my ($fasta_ref, undef, $header_ref) = @{read_fasta_file($fasta)};
my $n=0;
while (my ($id, $seq) = each %{$fasta_ref}) {

	# print out how many stuff we already computed
	$n++;
	if ($n % 100 == 0) {
		print STDERR '.';
		if ($n % 1000 == 0) {
			print STDERR $n / 1000, 'K';
		}
	}

	my $affinity = $header_ref->{$id};
	say join(' ', 't', $id, $affinity);
	my $graph_id = 0;

	# create sequence edges and vertices
	my @seq = split(//, $seq);
	foreach my $nt (@seq) {
		say join(' ', 'v', $graph_id, $nt, '1');
		if ($graph_id > 0) {
			say join(' ', 'e', $graph_id-1, $graph_id, 'b'); # b like backbone
		}
		$graph_id++;
	}

	# skip accessibility part if requested
	next if ($nostruct);

	# compute accessibilities on the fly
	my @accs = @{computeAccessibilities($seq)};

#	# get sequence accessibilities from file [DEPRECATED}
#	my @accs = (readAccessibilities($seq,$path,'', 1),
#		readAccessibilities($seq,$path,'_E', 1),
#		readAccessibilities($seq,$path,'_H', 1),
#		readAccessibilities($seq,$path,'_I', 1),
#		readAccessibilities($seq,$path,'_M', 1));

	# for each seq position sort context acc labels
	# and generate graph
	for (my $pos = 0; $pos < length($seq); $pos++) {

		# at first create the structure context subgraph
		my %accs;
		if ($nocontext) {
		    %accs = (
		        'U' =>$accs[0][$pos]);
		} else {
		    %accs = (
			    'U' => 1, # ensure that label for unpaired is always on top
			    'E' => $accs[1][$pos],
			    'H' => $accs[2][$pos],
			    'I' => $accs[3][$pos],
			    'M' => $accs[4][$pos]);
		}
		
		foreach my $key (keys %accs) {
			if ($accs{$key} <= $cutoff) {
				delete $accs{$key};
			}
		}

		my $num_vertices = keys %accs;
		my $unpaired_prob = $accs[0][$pos];
		if ($num_vertices > 1 or $unpaired_prob > $cutoff) {
			# if unpaired values above cutoff exist do the whole shebang
			
			my @keys = sort { $accs{$b} <=> $accs{$a} } keys %accs;
			my @node_ids = ($pos, $graph_id..($graph_id+$num_vertices-1));
			# save id of unpaired vertice for later
			my $unpaired_id = $graph_id;
			$graph_id += $num_vertices;

			while (@keys) {
				say join(' ', 'v', $node_ids[1], $keys[0], '0');
				# don't link the first vertice
				if (@keys != $num_vertices) {
					say join(' ', 'e', $node_ids[0], $node_ids[1], 'a'); # a like accessibility
				}
				shift @keys;
				shift @node_ids;
			}

			# and connect the subgraphs depending on the relation
			# of the paired and unpaired probabilities
			my $paired_prob = 1-$unpaired_prob;
			if ($paired_prob > $cutoff) {
				# do it the normal way
				
				# create the paired probability vertice
				my $paired_id = $graph_id;
				say join(' ', 'v', $paired_id, 'P', '0');
				$graph_id++;

				if ($unpaired_prob > $paired_prob) {
					say join(' ', 'e', $pos, $unpaired_id, 'a'); # a like accessibility
					say join(' ', 'e', $unpaired_id, $paired_id, 'a'); # a like accessibility
				}
				else {
					say join(' ', 'e', $pos, $paired_id, 'a'); # a like accessibility
					say join(' ', 'e', $paired_id, $unpaired_id, 'a'); # a like accessibility
				}
			} else {
				# leave out paired probability
				say join(' ', 'e', $pos, $unpaired_id, 'a'); # a like accessibility
			}
		} else {
			# else just print out the paired probability
			
			# create the paired probability vertice
			my $paired_id = $graph_id;
			say join(' ', 'v', $paired_id, 'P', '0');
			$graph_id++;
			
			say join(' ', 'e', $pos, $paired_id, 'a'); # a like accessibility
		}
	}
}

