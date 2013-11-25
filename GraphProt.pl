#!/usr/bin/perl
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;
use File::Copy;

=head1 NAME

GraphProt.pl

=head1 SYNOPSIS

GraphProt.pl -mode {regression,classification} -action {ls,train,test,cv,ntmargins,motif}

Options:

    -mode        'regression' or 'classification'
                     default: classification
    -action      what should GraphProt do?
                     ls: optimize parameters
                     cv: run a crossvalidation
                     train: train a model
                     predict: predict margins given a model
                     predict_nt: predict nucleotide-wise margins given a model
                     motif: create sequence and structure motifs given a model
    -onlyseq     use GraphProt sequence models
    -fasta       fasta file containing binding sites
    -affinities  list of affinities
                     one value per line, same order as binding sites (fasta)
    -negfasta    fasta file containing negative class sequences
    -abstraction RNAshapes abstraction level [RNA structure graphs]
                     default: 3
    -R           GraphProt radius
                     default: 1
    -D           GraphProt distance
                     default: 4
    -bitsize     GraphProt bitsize used for feature encoding
                     default: 14
    -c           SVR parameter c       [regression]
                     default: 1
    -epsilon     SVR parameter epsilon [regression]
                     default: 0.1
    -lambda      SGD parameter lambda  [classification]
                     default: 10e-4
    -epochs      SGD parameter epochs  [classification]
                     default: 10
    -debug       enable debug output
    -help        brief help message
    -man         full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# create temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
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
# defaults
###############################################################################
my $def_abstraction = 3;
my $def_R = 1;
my $def_D = 4;
my $def_bitsize = 14;
my $def_epochs = 10;
my $def_lambda = 10e-4;
my $def_epsilon = 0.1;
my $def_c = 1;

###############################################################################
# check program paths
###############################################################################

# TODO check RNAshapes
# TODO check EDeN
# TODO check make

###############################################################################
# parse command line options
###############################################################################
my $mode;
my $action;
my $onlyseq;
my $fasta;
my $negfasta;
my $affys;
my $R;
my $D;
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
                            "R=i" => \$R,
                            "D=i" => \$D,
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

sub check_param_fasta {
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

sub check_param_negfasta {
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

sub check_param_affys {
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

sub check_param_R {
#    if (not defined $R) {
#        say STDERR "missing parameter: please specify the radius (R)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $R or $R = $def_R;
    if ($R < 0) {
        pod2usage("error: please specify a positive radius (R)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_D {
#    if (not defined $D) {
#        say STDERR "missing parameter: please specify the distance (D)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $D or $D = $def_D;
    if ($D < 0) {
        pod2usage("error: please specify a positive distance (D)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_bitsize {
#    if (not defined $bitsize) {
#        say STDERR "missing parameter: please specify the bitsize (bitsize)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $bitsize or $bitsize = $def_bitsize;
    if ($bitsize < 8) {
        pod2usage("error: please specify a positive bitsize larger than 8 (bitsize)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_abstraction {
#    if (not defined $abstraction) {
#        say STDERR "missing parameter: please specify the RNAshapes abstraction level (abstraction)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $abstraction or $abstraction = $def_abstraction;
    if ($abstraction < 1 or $abstraction > 5) {
        pod2usage("error: please specify an RNAshapes abstraction level between 1 and 5 (abstraction)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_c {
#    if (not defined $c) {
#        say STDERR "missing parameter: please specify SVR parameter c (c)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $c or $c = $def_c;
    if ($c <= 0) {
        pod2usage("error: please specify the SVR parameter c (c)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_epsilon {
#    if (not defined $epsilon) {
#        say STDERR "missing parameter: please specify SVR parameter epsilon (epsilon)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $epsilon or $epsilon = $def_epsilon;
    if ($c <= 0) {
        pod2usage("error: please specify a positive epsilon (epsilon)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_epochs {
#    if (not defined $epochs) {
#        say STDERR "missing parameter: please specify SGD parameter epochs (epochs)";
#        $quit_with_help = 1;
#        return 0;
#    }
    defined $epochs or $epochs = $def_epochs;
    if ($epochs <= 0) {
        pod2usage("error: please specify a value larger 0 (epochs)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_param_lambda {
    defined $lambda or $lambda = $def_lambda;
#    if (not defined $lambda) {
#        say STDERR "missing parameter: please specify SGD parameter lambda (lambda)";
#        $quit_with_help = 1;
#        return 0;
#    }
    if ($lambda <= 0) {
        pod2usage("error: please specify a positive lambda (lambda)");
        $quit_with_help = 1;
        return 0;
    }
}

sub check_params_regression {
    # only check when creating a structure model
    defined $onlyseq or check_param_R();
    defined $onlyseq or check_param_D;
    defined $onlyseq or check_param_abstraction;
    check_param_bitsize;
    check_param_c;
    check_param_epsilon;
}

sub check_params_classification {
    # only check when creating a structure model
    defined $onlyseq or check_param_R;
    defined $onlyseq or check_param_D;
    defined $onlyseq or check_param_abstraction;
    check_param_bitsize;
    check_param_epochs;
    check_param_lambda;
}

defined $action or pod2usage("please specify the GraphProt action\n");

if ($mode eq 'regression') {
    if ($action eq 'ls') {
        check_param_fasta;
        check_param_affys;
    } elsif ($action eq 'cv') {
        check_param_fasta;
        check_param_affys;
        check_params_regression;
    } elsif ($action eq 'train') {
        check_param_fasta;
        check_param_affys;
        check_params_regression;
    } elsif ($action eq 'predict') {
        check_param_fasta;
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
        check_param_fasta;
        check_param_negfasta;
    } elsif ($action eq 'cv') {
        check_param_fasta;
        check_param_negfasta;
        check_params_classification;
    } elsif ($action eq 'train') {
        check_param_fasta;
        check_param_negfasta;
        check_params_classification;
    } elsif ($action eq 'predict') {
        check_param_fasta;
        check_param_negfasta;
        check_params_classification;
    } elsif ($action eq 'predict_nt') {
        check_param_fasta;
        check_params_classification;    
    } elsif ($action eq 'motif') {
        check_param_fasta;
        check_params_classification;    
    } else {
        pod2usage("error: unknown action '$action'\n");
    }
} else {
    pod2usage("error: unknown mode '$mode'\n");
}

$quit_with_help and pod2usage();


###############################################################################
# main
###############################################################################

# set up temporary directory
my $tmp_template = 'GraphProt-XXXXXX';
my $tmp_prefix = './';
my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 0);

# write parameters
if ($action ne "ls") {
    open PARAMETERS, ">", "$tmpdir/ID.param";
    defined $R and say PARAMETERS "R $R";
    defined $D and say PARAMETERS "D $D"; 
    defined $c and say PARAMETERS "c $c"; 
    defined $epsilon and say PARAMETERS "e $epsilon"; 
    defined $epochs and say PARAMETERS "EPOCHS $epochs"; 
    defined $lambda and say PARAMETERS "LAMBDA $lambda"; 
    defined $abstraction and say PARAMETERS "ABSTRACTION $abstraction"; 
    defined $bitsize and say PARAMETERS "b $bitsize"; 
    say PARAMETERS "CUE nil";
    say PARAMETERS "STACK nil";
    say PARAMETERS "VIEWPOINT --vp";
    close PARAMETERS;
}

# copy input files
# TODO

# collect make call
my $scriptdir = dirname($0);
my $makecall = "make -C $scriptdir ";

# use sequence graphs
if (defined $onlyseq) {
    $makecall .= " -e GRAPH_TYPE=SEQUENCE"
} else {
# use structure graphs
    $makecall .= " -e GRAPH_TYPE=CONTEXTSHREP"
}

if ($mode eq 'regression') {
    if ($action eq 'ls') {
        # add parameters
        $makecall .= " -e SVM=SVR ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'cv') {
        # add parameters
        $makecall .= " -e SVM=SVR ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'train') {
        # add parameters
        $makecall .= " -e SVM=SVR ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'predict') {
        # add parameters
        $makecall .= " -e SVM=SVR ";
        # add targets
        $makecall .= " ";
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
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'cv') {
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'train') {
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'predict') {
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'predict_nt') {
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } elsif ($action eq 'motif') {
        # add parameters
        $makecall .= " -e SVM=SGD ";
        # add targets
        $makecall .= " ";
    } else {
        pod2usage("error: unknown action '$action'\n");
    }
} else {
    pod2usage("error: unknown mode '$mode'\n");
}

# execute make call
# get output files
# clean up
