# user parameters

# written out with end results
EXPERIMENT_ID:=default

# proteins to use
#PROTEINS:=Fusip HuR PTB RBM4 SF2 SLM2 U1A VTS1 YB1
PROTEINS:=test

# sets for training: FULL/WEAK/ALL
TRAINING_SETS=FULL

# enable/disable line search: YES/NO
DO_LINESEARCH=YES

# set graph type:
# ONLYSEQ: plain sequence
# STRUCTACC: sequence annotated with structural context
# NORMSHREP: plain shrep graphs
# PROBSHREP: plain shrep graphs, features weighted by shape probability
# CONTEXTSHREP: shrep graphs annotated with structural context
GRAPH_TYPE=ONLYSEQ

# main
SHELL:=/bin/bash
.DELETE_ON_ERROR:

# don't delete intermediate files
.SECONDARY:

# parameters:
LINESEARCH_INPUT_SIZE:=5000
SVR_CACHE:=10000
CV_FOLD:=5

# paths
ROOT:=~/projects/RBPaffinity
FA_DIR:=~$(ROOT)/data/fasta
THR_DIR:=$(ROOT)/data/thresholds/

# binaries
SVRTRAIN:=~/src/libsvm-3.0/svm-train -s 3 -t 0 -m $(SVR_CACHE)
SVRPREDICT:=~/src/libsvm-3.0/svm-predict
PERF:=~/src/stat/perf
LINESEARCH:=./lineSearch.pl
COMBINEFEATURES:=./bin/combineFeatures.pl
SHUF:=~/src/coreutils-8.15/src/shuf
FASTAPL:=/usr/local/user/RNAtools/fastapl
FASTA2GSPAN:=/usr/local/user/RNAtools/fasta2shrep_gspan.pl
NSPDK:=~/projects/RBPaffinity/bin/NSPDK

# targets
FULL_BASENAMES:=$(patsubst %,%_data_full_A,$(PROTEINS)) \
			$(patsubst %,%_data_full_B,$(PROTEINS))
BRUIJN_BASENAMES:=$(patsubst %,%_data_bruijn_A,$(PROTEINS)) \
			$(patsubst %,%_data_bruijn_B,$(PROTEINS))
ifeq ($(TRAINING_SETS),FULL)
BASENAMES:=$(FULL_BASENAMES)
else
ifeq ($(TRAINING_SETS),WEAK)
BASENAMES:=$(BRUIJN_BASENAMES)
else
BASENAMES:=$(FULL_BASENAMES) $(BRUIJN_BASENAMES)
endif
endif

PERF_FILES:=$(patsubst %,%.perf,$(BASENAMES))
PARAM_FILES:=$(patsubst %,%.param,$(BASENAMES))
CV_FILES:=$(patsubst %,%.cv,$(BASENAMES))
CSTAT_FILES:=$(patsubst %,%.cstats,$(FULL_BASENAMES))

# evaluations specific to graph type
################################################################################
ifeq ($(GRAPH_TYPE),ONLYSEQ)
# line search parameters
LSPAR:=./ls.onlyseq.parameters

%.gspan : %.fa
	/usr/local/perl/bin/perl $(ROOT)/bin/create_accgraph/createExtendedGraph.pl \
	--nostruct -fa $< > $@

# feature creation for this type of graph has to set an additional parameter
# in order to center only on nucleotides and not on annotation -> -T nspdkvp
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	# create features
	ln -sf $< $* # remove suffix to have shorter filenames
	$(NSPDK) -fg $* -of -R $(RADIUS) -D $(DISTANCE) -b $(BITSIZE) -T nspdkvp -gt $(DIRECTED)
	-rm -f $* $@_bin # clean up after feature creation
	# add affinities to features
	mv $@ $@.tmp
	cat $@.tmp | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -rf $@.tmp # clean up affinityless feature file

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$NF}' > $@

else
################################################################################
ifeq ($(GRAPH_TYPE),STRUCTACC)
# line search parameters
LSPAR:=./ls.structacc.parameters

%.gspan : %.fa
	time /usr/local/perl/bin/perl ./bin/create_accgraph/createExtendedGraph.pl \
	-fa $< > $@

# feature creation for this type of graph has to set an additional parameter
# in order to center only on nucleotides and not on annotation -> -T nspdkvp
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	# create features
	ln -sf $< $* # remove suffix to have shorter filenames
	$(NSPDK) -fg $* -of -R $(RADIUS) -D $(DISTANCE) -b $(BITSIZE) -T nspdkvp -gt $(DIRECTED)
	-rm -f $* $@_bin # clean up after feature creation
	# add affinities to features
	mv $@ $@.tmp
	cat $@.tmp | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -rf $@.tmp # clean up affinityless feature file

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$NF}' > $@

else
################################################################################
ifeq ($(GRAPH_TYPE),NORMSHREP)
# line search parameters
LSPAR:=./ls.shrep.parameters

# else, create gspan from fasta
%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	time $(FASTA2GSPAN) $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -fasta $< > $@

# feature creation for this type of graph has to set an additional parameter
# in order to center only on nucleotides and not on annotation -> -T nspdkvp
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : bitsize=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	ln -sf $< $* # remove suffix to have shorter filenames
	$(ROOT)/bin/NSPDK -fg $* -of -R $(RADIUS) -D $(DISTANCE) -b $(bitsize) -T nspdkvp -gt DIRECTED
	-rm -f $* $@_bin # clean up after feature creation
	mv $@ $@.tmp
	cat $@.tmp | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -rf $@.tmp # clean up affinityless feature file

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

else
################################################################################
ifeq ($(GRAPH_TYPE),PROBSHREP)
# line search parameters
LSPAR:=./ls.shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	time $(FASTA2GSPAN) $(STACK) $(CUE) -stdout -q -Tp 0.05 -t $(ABSTRACTION) -M 5 -fasta $< > $@

# feature creation for this type of graph has to set an additional parameter
# in order to center only on nucleotides and not on annotation -> -T nspdkvp
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : bitsize=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	# create features
	# remove t and w, convert s to t
	cat $< | grep -v -e '^t' -e '^w' | sed 's/^s/t/' > $*_singleshreps
	$(ROOT)/bin/NSPDK -fg $*_singleshreps -of -R $(RADIUS) -D $(DISTANCE) -b $(bitsize) -gt $(DIRECTED)
	-rm -f $*_singleshreps $*_singleshreps.feature_bin # clean up after feature creation
	# write out probabilities
	cat $< | grep '^s' | perl -ne '($$prob) = /SHAPEPROB ([0-9.]+)/; print $$prob, "\n"' > $*_probs
	# write out shrep membership
	cat $< | awk '/^t/{i++}/^s/{print i}' > $*_groups
	# go go perl
	/usr/local/perl/bin/perl $(COMBINEFEATURES) $*_singleshreps.feature $*_probs $*_groups > $*
	-rm -f $*_singleshreps.feature $*_probs $*_groups
	# add affinities to features
	cat $* | grep -v \"^\$\" | paste -d' ' $*.affy - > $@

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

else
################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=./ls.shrep.parameters
endif
endif
endif
endif
endif

.PHONY: all ls cv classstats

all: $(PERF_FILES) results_aucpr.csv

ls : $(PARAM_FILES)

cv : $(CV_FILES)

classstats : summary.cstats $(CSTAT_FILES)

%.cv : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv : %.feature %.param
	$(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# if available, create gspan from precomputed files
%.gspan : %.gspan.gz
	zcat $< > $@

# link parameter file for simpler handling
%.pred.param : %.param
	ln -sf $< $@

ifeq ($(DO_LINESEARCH),NO)
# just use defaults instead of doing line search
%.param : $(LSPAR)
	cut -f 1,2 -d' ' < $< > $@
else
# do parameter optimization by line search
%.param : %.ls.fa $(LSPAR)
	$(LINESEARCH) -fa $< -param $(LSPAR) -mf Makefile -of $@ 2> >(tee $@.log >&2)
endif

# subset original fasta
%.ls.fa : %.fa
	cat $< | \
	$(FASTAPL) -e 'print ">", $$head, "\t", $$seq, "\n"' | \
	$(SHUF) -n $(LINESEARCH_INPUT_SIZE) | \
	/usr/local/perl/bin/perl -ane \
	'$$seq = pop @F; $$head = join(" ", @F); print $$head, "\n", $$seq, "\n";' > \
	$@

%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature %.param
	$(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

%.svrout : %.model %.pred.feature
	$(SVRPREDICT) $*.pred.feature $< $@

%.pred : %.svrout %.pred.affy
	# combine affinities and predictions
	paste $*.pred.affy $< > $@

%.perf : BASENAME=$(firstword $(subst _, ,$<))
%.perf : $(shell echo $(BASENAME))
%.perf : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.perf : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.perf : %.pred
	# select by threshold
	( cat $< | awk '$$1 > $(HT)' | cut -f 2 | awk '{print 1 "\t" $$1 }'; \
	cat $< | awk '$$1 < $(LT)' | cut -f 2 | awk '{print 0 "\t" $$1}' ) > $@.threshold
	# compute performance measures
	$(PERF) < $@.threshold > $@
	rm -rf $@.threshold*

results_aucpr.csv : $(PERF_FILES)
	grep ROC $(PERF_FILES) | tr ':' ' ' | \
	awk '{print $$1, "$(EXPERIMENT_ID)", $$NF}' | sort > roc.tmp
	grep APR $(PERF_FILES) | tr ':' ' ' | \
	awk '{print $$1, "$(EXPERIMENT_ID)", $$NF}' | cut -f 1,3 -d' ' | \
	sort > aucpr.tmp
	join roc.tmp aucpr.tmp > $@
	rm -rf roc.tmp aucpr.tmp

%.cstats : BASENAME=$(firstword $(subst _, ,$<))
%.cstats : TYPE=$(word 3,$(subst _, ,$<))
%.cstats : SET=$(word 4,$(subst ., ,$(subst _, ,$<)))
%.cstats : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.cstats : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.cstats : HN=$(shell cat $< | grep '^>' | awk '$$NF > $(HT)' | wc -l)
%.cstats : LN=$(shell cat $< | grep '^>' | awk '$$NF < $(LT)' | wc -l)
%.cstats : %.fa
	perl -e 'print join("\t", "$(BASENAME)", "$(SET)", "$(LT)", "$(LN)", "$(HT)", "$(HN)"),"\n"' > $@

summary.cstats : $(CSTAT_FILES)
	( perl -e 'print join("\t", "protein", "set", "negative threshold", "negative instances", "positive threshold", "positive instances"),"\n"'; \
	cat $^ | sort -k1,2 ) > $@

# keep predictions and results
clean:
	-rm -rf $(MODELS) log *.gspan *.threshold* *.model *.feature *.affy *.svrout *.ls.fa *.log

# delete all files
distclean: clean
	-rm -rf *.param *.fa *.perf *.pred
