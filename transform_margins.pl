#!/usr/local/perl/bin/perl
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

transform_margins.pl -W WINDOW_SIZE

reads nucleotide-wise margins from stdin or file,
computes min, max, mean, median, sum per window of size WINDOW_SIZE

input format (tab separated): sequence id, sequence position, margin
output is input plus: min, max, mean, median, sum

=head1 SYNOPSIS

Options:

    -W          window size
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# create temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
#my $tmp_template = 'template-XXXXXX';
#my $tmp_prefix = '/var/tmp/';
#my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);
#$SIG{'INT'} = 'end_handler';
#$SIG{'TERM'} = 'end_handler';
#$SIG{'ABRT'} = 'end_handler';
#$SIG{'USR1'} = 'end_handler';
#$SIG{'USR2'} = 'end_handler';
#sub end_handler {
#	print STDERR "signal '", $_[0], "' caught, cleaning up temporary files\n";
#	# change into home directory. deletion of the temporary directory will
#	# fail if it is the current working directory
#	chdir();
#	File::Temp::cleanup();
#	die();
#}

###############################################################################
# parse command line options
###############################################################################
my $winsize;
my $help;
my $man;
my $debug;
my $result = GetOptions (	"help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug,
                            "W=i"   => \$winsize);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

(defined $winsize) or pod2usage("error: parameter W mandatory");
($winsize>0) or pod2usage("error: window size must be larger 0");

###############################################################################
# functions
###############################################################################

sub summarize {
  return 'dummyvalue';
}

sub cmp_wins {
  # array of lines to summarize
  my ($linestack_aref) = @_;
  my @lines = @{$linestack_aref};
  # array of lines containing summarized values
  my @summaries;
  # queue of margins to summarize (the window)
  my @winqueue;

  # handle left border and full windows
  foreach my $line (@lines) {
    # push margin into queue
    my (undef, undef, $margin) = split("\t", $line);
    push @winqueue, $margin;
    # if window size reached, remove oldest item from queue
    shift @winqueue if (@winqueue > $winsize);
    # compute summaries if at least ceil(W/2) items in queue
    if (@winqueue >= ceil($winsize/2)) {
        # summarize window values
        push @summaries, summarize(\@winqueue);
    }
  }

  # handle remaining windows on right border
  while (shift @winqueue) {
    push @summaries, summarize(\@winqueue) if (@winqueue >= floor($winsize/2)+1)
  }

  # print output
  my $nlines = @lines;
  my $nsummaries = @summaries;
  say STDERR "found $nsummaries summaries for $nlines lines";
  foreach my $line (@lines) {
    my $summary = shift @summaries;
    say join("\t", $line, $summary);
  }

}

###############################################################################
# main
###############################################################################

# parse input; data for each sequence is forwarded to the windowing function
my $current_seqid = 0; # fist sequence id is expected to be 0
my @linestack;
while (my $line = <>) {
  chomp $line;
  my ($seqid) = split("\t", $line);
  if ($seqid == $current_seqid) {
    # if sequence id not changed, just record line
    push @linestack, $line;
  } else {
    # push lines to windowing function
    cmp_wins(\@linestack);
    # start new linestack using current line
    @linestack = ($line);
    # update current_seqid
    $current_seqid = $seqid;
  }
}
# push lines to windowing function
cmp_wins(\@linestack)
