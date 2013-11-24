#!/usr/bin/perl
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;

=head1 NAME

GraphProt.pl

=head1 SYNOPSIS

GraphProt.pl -mode {regression,classification} -action {ls,train,test,cv,ntmargins,motif}

Options:

    -mode       'regression' or 'classification'
                    default: classification
    -action     what should GraphProt do?
                    ls: optimize parameters
                    cv: run a crossvalidation
                    train: train a model
                    predict: predict margins given a model
                    predict_nt: predict nucleotide-wise margins given a model
                    motif: create sequence and structure motifs given a model
    -onlyseq    use GraphProt sequence models
    -fasta      fasta file containing binding sites
    -affinities list of affinities
                    one value per line, same order as binding sites (fasta)
    -negfasta   fasta file containing negative class sequences
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
my $tmp_template = 'template-XXXXXX';
my $tmp_prefix = '/var/tmp/';
my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);
$SIG{'INT'} = 'end_handler';
$SIG{'TERM'} = 'end_handler';
$SIG{'ABRT'} = 'end_handler';
$SIG{'USR1'} = 'end_handler';
$SIG{'USR2'} = 'end_handler';
sub end_handler {
	print STDERR "signal '", $_[0], "' caught, cleaning up temporary files\n";
	# change into home directory. deletion of the temporary directory will
	# fail if it is the current working directory
	chdir();
	File::Temp::cleanup();
	die();
}

###############################################################################
# parse command line options
###############################################################################
my $mode;
my $action;
my $onlyseq;
my $fasta;
my $negfasta;
my $affys;
my $radius;
my $distance;
my $c;
my $epsilon;
my $epochs;
my $lambda;
my $abstraction;
my $bitsize;
my $help;
my $man;
my $debug;
my $result = GetOptions (	"mode=s" => \$mode,
                            "action=s" => \$action,
                            "onlyseq" => \$onlyseq,
                            "fasta=s" => \$fasta,
                            "negfasta=s" => \$negfasta,
                            "affinities=s" => \$affys,
                            "R=i" => \$radius,
                            "D=i" => \$distance,
                            "c=f" => \$c,
                            "epsilon=f" => \$epsilon,
                            "epochs=i" => \$epochs,
                            "lambda=f" => \$lambda,
                            "abstraction=i" => \$abstraction,
                            "bitsize=i" => \$bitsize,
                            "help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

###############################################################################
# check parameters
###############################################################################

(defined $mode) or $mode='classification';
(defined $action) or pod2usage("please specify -action\n");

my $quit_with_help = 0;

sub check_fasta {
    if (not defined $fasta) {
        say STDERR "missing parameter: please specify a set of binding sites (fasta)";
        $quit_with_help = 1;
        return 0;
    }
    if (not -f $fasta) {
        pod2usage("error: can't find file '$fasta'");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_negfasta {
    if (not defined $negfasta) {
        say STDERR "missing parameter: please specify a set of unbound sites";
        $quit_with_help = 1;
        return 0;
    }
    if (not -f $negfasta) {
        pod2usage("error: can't find file '$negfasta'");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_affys {
    if (not defined $affys) {
        say STDERR "missing parameter: please specify the list of affinities (affinities)";
        $quit_with_help = 1;
        return 0;
    }
    if (not -f $affys) {
        pod2usage("error: can't find file '$affys'");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_params_regression {
    $onlyseq or check_R
    $onlyseq or check_D
    $onlyseq or check_abstraction
    check_bitsize
    check_c
    check_epsilon
}

sub check_params_classification {
    $onlyseq or check_R
    $onlyseq or check_D
    $onlyseq or check_abstraction
    check_bitsize
    check_epochs
    check_epsilon
}

defined $action or pod2usage("please specify the GraphProt action\n");

if ($mode eq 'regression') {
    if ($action eq 'ls') {
        check_fasta;
        check_affys;
    } elsif ($action eq 'cv') {
        check_fasta;
        check_affys;
        check_params_regression;
    } elsif ($action eq 'train') {
        check_fasta;
        check_affys;
        check_params_regression;
    } elsif ($action eq 'predict') {
        check_fasta;
        check_params_regression;
    } elsif ($action eq 'predict_nt') {
        say STDERR "sorry, invalid action in regression setting";
        exit 2;
    } elsif ($action eq 'motif') {
        say STDERR "sorry, invalid action in regression setting";
        exit 2;
    } else {
        pod2usage("error: unknown action '$action'\n");
    }
} elsif ($mode eq 'classification') {
    if ($action eq 'ls') {
        check_fasta;
        check_negfasta;
    } elsif ($action eq 'cv') {
        check_fasta;
        check_negfasta;
        check_params_classification;
    } elsif ($action eq 'train') {
        check_fasta;
        check_negfasta;
        check_params_classification;
    } elsif ($action eq 'predict') {
        check_fasta;
        check_negfasta;
        check_params_classification;
    } elsif ($action eq 'predict_nt') {
        check_fasta;
        check_params_classification;    
    } elsif ($action eq 'motif') {
        check_fasta;
        check_params_classification;    
    } else {
        pod2usage("error: unknown action '$action'\n");
    }
} else {
    pod2usage("error: unknown mode '$mode'\n");
}

$quit_with_help and pod2usage();

###############################################################################
# check program paths
###############################################################################


###############################################################################
# main
###############################################################################

