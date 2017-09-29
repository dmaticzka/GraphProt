#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;

=head1 NAME

=head1 SYNOPSIS

toDot.pl < [gpsan file]

directly plots graph to graph_id.png

Options:

	-plotno     plot the nth graph in input file (default 1)
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
# getColor
# in: label
# out: node attributes
###############################################################################
my %ncol = ("P" => ", color=red, style=filled",
			"A" => ", color=red",
			"T" => ", color=blue",
			"C" => ", color=green",
			"G" => ", color=orange");
sub getAttr {
    my ($label) = @_;
    my $color = $ncol{$label};
    (defined $color) or $color='black';
    return $color;
}

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $plotno;
my $result = GetOptions (	"help"	=> \$help,
							"man"	=> \$man,
							"plotno=i" => \$plotno);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);
(defined $plotno) or $plotno = 1;

###############################################################################
# main
###############################################################################
my $i;
while (<>) {
	if (/^t/) {
		$i++;
		if ($plotno == $i) {
			my (undef, $id, $affy) = split();
            ($id !~ /#/) or $id = $i;
			open DOT, ">tmp.dot";
#            open DOT, ">&STDOUT";
			say DOT "digraph $id {";
#			say DOT "affinity [label=$affy];";
		}
	} elsif (/^v/ and ($plotno == $i)) {
		my (undef, $id, $label) = split();
		$label =~ s/#//;
		$label =~ s/\^//;
		my $attr = getAttr($label);
		say DOT "$id [label=\"$label $attr\"];";
	} elsif (/^e/ and ($plotno == $i)) {
		my (undef, $first, $second, $label) = split();
		say DOT "$first -> $second;";
	}
	if ($i > $plotno) {
		say DOT '}';
		close DOT;
		exit;
	}
}
if (not $i > $plotno) {
		say DOT '}';
		close DOT;
		exit;
}
