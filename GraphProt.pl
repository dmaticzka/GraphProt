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

=head1 SYNOPSIS

Options:

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
my $help;
my $man;
my $debug;
my $result = GetOptions (	"help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################

