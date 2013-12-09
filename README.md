# GraphProt Manual #

This software package contains the GraphProt framework as published in
"GraphProt: modeling binding preferences of
RNA-binding proteins".

## Installation ##

GraphProt contains a precompiled version of "EDeN", the SVM package used for
feature creation and classification. This binary should run on most Linux-based
systems. In case it does not run on your system, please call
"bash ./recompile_EDeN.sh" from the GraphProt main directory.

GraphProt uses various opensource software packages. Please make sure that the
follwing programs are installed and accessible via the PATH environment variable
(i.e. you should be able to call the programs by just issuing the command).

* RNAshapes is used for GraphProt secondary structure predictions (recommended version: 2.1.6, [http://bibiserv.techfak.uni-bielefeld.de/rnashapes/](http://bibiserv.techfak.uni-bielefeld.de/rnashapes/))
* perf is used to calculate prediction performance ([http://osmot.cs.cornell.edu/kddcup/software.html](http://osmot.cs.cornell.edu/kddcup/software.html))
* libsvm is used for support vector regressions ([http://www.csie.ntu.edu.tw/~cjlin/libsvm/](http://www.csie.ntu.edu.tw/~cjlin/libsvm/))
* GNU make is used as the pipeline backend ([http://www.gnu.org/software/make/](http://www.gnu.org/software/make/))
* R is used to process nucleotide-wise margins for motif creation ([www.r-project.org/](www.r-project.org/))
* The R plyr package is required for calculating motifs ([http://plyr.had.co.nz/](http://plyr.had.co.nz/)) and can be installed from within R by issuing the command "install.packages('plyr')".

GraphProt will scan for these programs and notify you if something seems amiss.
GraphProt contains a copy of fastapl ([http://seq.cbrc.jp/fastapl/index.html.en](http://seq.cbrc.jp/fastapl/index.html.en)).

## Usage ##

GraphProt analyses are started by calling "GrapProt.pl". If no options are given,
GraphProt.pl will display a help message summarizing all available options.
The default mode is to run analyses in classification setting,
switch to regression setting using the parameter -mode regression.
In general, GraphProt analyses are run by issuing different actions, e.g.

  GraphProt.pl --action train -fasta train_positives.fa -negfasta train_negatives.fa

GraphProt supports input sequences in fasta format. The **viewpoint** mechanism
sets viewpoints to all nucleotides in **uppercase** letters, nucleotides in
**lowercase** letters are only used for RNA structure predictions.

GraphProt parameters abstraction, R, D, bitsize, c, epsilon, epochs and lambda
are set to default values. For best results, optimized parameters should be
obtained with the ls parameter optimization setting.

Input files in classification setting are specified with parameters "-fasta"
(binding sites) and "-negfasta" (unbound sites). For regressions, input sequences
are specified with "-fasta" and sequence scores with "-affinities". For each
sequence, the affinity file should contain one value per line.

Output filenames can be specified via a prefix (-prefix); if no prefix is given,
the default is "GraphProt".

### Available Actions ###

#### ls - Parameter Optimization ####

Determines optimized parameters. Parameters are printed to screen and written
to file "GraphProt.param".

#### cv - Crossvalidation ####

Runs a 10-fold crossvalidation. Crossvalidation results are written to file
"GraphProt.cv_results".

#### train - Model Training ####

Trains a GraphProt model. The model is written to file "GraphProt.model".

#### predict - Predictions based on a model ####

Predict SVM margins for whole sequences. Margins are written to file
"GraphProt.predictions". Each line of this file gives the margin for one sequence in the second column,
in the same order as the fasta file. In classification setting the first column contains the class,
in regression setting the first column contains, if specified, the affinities, otherwise
1.

#### predict_nt - Nucleotide-wise predictions based on a model ####

Predict nucleotide-wise margins for sequences. Nucleotide-wise margins are written
to file "GraphProt.nt_margins" and contain three columns:

1. number of sequence
2. number of nucleotide
3. margin of this nucleotide

Please note that nucleotide wise margins are only supported for sequences with
at most 150 nt.

#### motif - Create RNA sequence and structure motifs ####

Create RNA sequence and structure motifs. High-scoring 12-mers are written to
files "GraphProt.sequence_motif" and "GraphProt.structure_motif" and can be
used to create sequence logos (for example using WebLogo 
[http://weblogo.threeplusone.com/](http://weblogo.threeplusone.com/)).

## Advanced Usage ##

In addition to the integrated usage via GraphProt.pl, individual tasks such as
creation of RNA structure graphs or calculation of features can be accomplished
using the following tools:

* fasta2shrep_gspan.pl: graph creation
* EDeN/EDeN: NSPD kernel and SGD support vector machine

Usage information for these tools can be obtained by specifying the "-h" option.

### RNA sequence and structure graphs ###

RNA sequence and structure graphs are created using fasta2shrep_gspan.pl. Structure graphs
are created using the following parameters. The user has to chose an appropriate
RNAshapes __ABSTRACTION_LEVEL__.

  fasta2shrep_gspan.pl --seq-graph-t --seq-graph-alph -abstr -stdout -M 3 -wins '150,' -shift '25' -fasta PTBv1.train.fa -t __ABSTRACTION_LEVEL__ | gzip > PTBv1.train.gspan.gz
  
RNA sequence graphs are created using the following parameters:

  fasta2shrep_gspan.pl --seq-graph-t -nostr -stdout -fasta PTBv1.train.fa | gzip > PTBv1.train.gspan.gz

### NSPD kernel and SGD support vector machine ###

For example, 10-fold crossvalidation using EDeN is done via:

  EDeN/EDeN -a CROSS_VALIDATION -c 10 -i PTBv1.train.gspan.gz -t PTBv1.train.class -g DIRECTED -b __BIT_SIZE__ -r __RADIUS__ -d __DISTANCE__ -e __EPOCHS__ -l __LAMBDA__

and setting the appropriate parameters for BIT_SIZE, RADIUS, DISTANCE, EPOCHS
and LAMBDA.
