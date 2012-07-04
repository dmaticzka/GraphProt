# user parameters
include PARAMETERS

SHELL:=/bin/bash
.DELETE_ON_ERROR:

# don't delete intermediate files
.SECONDARY:

# parameters:
LINESEARCH_INPUT_SIZE:=5000
SVR_CACHE:=10000
CV_FOLD:=5

# paths
PROJDIR:=~/projects/RBPaffinity
FA_DIR:=$(PROJDIR)/data/fasta
THR_DIR:=$(PROJDIR)/data/thresholds/
# expect binaries to reside in pwd, otherwise this variable must be overwritten
BINDIR:=$(shell pwd)

# binaries
PERL:=/usr/local/perl/bin/perl
SVRTRAIN:=~/src/libsvm-3.0/svm-train -s 3 -t 0 -m $(SVR_CACHE)
SVRPREDICT:=~/src/libsvm-3.0/svm-predict
PERF:=~/src/stat/perf
LINESEARCH:=$(PERL) $(BINDIR)/lineSearch.pl
COMBINEFEATURES:=$(PERL) $(BINDIR)/combineFeatures.pl
SHUF:=~/src/coreutils-8.15/src/shuf
FASTAPL:=$(PERL) /usr/local/user/RNAtools/fastapl
#FASTAPL:=~/repositories/RNAtools/fastapl
FASTA2GSPAN:=$(PERL) /usr/local/user/RNAtools/fasta2shrep_gspan.pl
#FASTA2GSPAN:=~/repositories/RNAtools/fasta2shrep_gspan.pl
SVMSGDNSPDK:=~/repositories/svmsgdnspdk_dir_dev/svmsgdnspdk
CREATE_EXTENDED_ACC_GRAPH:=$(PERL) $(BINDIR)/create_accgraph/createExtendedGraph.pl
MERGE_GSPAN:=$(PERL) $(BINDIR)/merge_gspan.pl

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
LSPAR:=./ls.structacc.parameters

%.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) --nostruct -fa $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$NF}' > $@

endif

################################################################################
ifeq ($(GRAPH_TYPE),STRUCTACC)
# line search parameters
LSPAR:=./ls.structacc.parameters

%.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) -fa $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$NF}' > $@

endif

################################################################################
ifeq ($(GRAPH_TYPE),SHREP)
# line search parameters
LSPAR:=./ls.shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -fasta $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

endif

################################################################################
ifeq ($(GRAPH_TYPE),PROBSHREP)
# line search parameters
LSPAR:=./ls.shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) -stdout -q -Tp 0.05 -t $(ABSTRACTION) -M 5 -fasta $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	# remove t and w, convert s to t
	cat $< | grep -v -e '^t' -e '^w' | sed 's/^s/t/' > $*_singleshreps
	# write out probabilities
	cat $< | grep '^s' | $(PERL) -ne '($$prob) = /SHAPEPROB ([0-9.]+)/; print $$prob, "\n"' > $*_probs
	# write out shrep membership
	cat $< | awk '/^t/{i++}/^s/{print i}' > $*_groups
	# compute features
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $*_singleshreps -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*_singleshreps
	# compute probability-weighted features
	$(COMBINEFEATURES) R$(RADIUS)D$(DISTANCE)$*_singleshrepsoutput.vec $*_probs $*_groups > $*
	# add affinities to features
	cat $* | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*_singleshrepsoutput.vec

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

endif

################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=./ls.shrep_context.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) -abstr -stdout -t $(ABSTRACTION) -M 5 -fasta $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : bitsize=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : NHF=$(shell grep '^NHF ' $*.param | cut -f 2 -d' ')
%.feature : RR=$(shell grep '^RR ' $*.param | cut -f 2 -d' ')
%.feature : RD=$(shell grep '^RD ' $*.param | cut -f 2 -d' ')
%.feature : RW=$(shell grep '^RW ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -kt ABSTRACT -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -anhf $(NHF) -rR $(RR) -rD $(RD) -rW $(RW) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

%.affy : %.gspan
	# extract affinities from gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

endif

################################################################################
ifeq ($(GRAPH_TYPE),MEGA)
# line search parameters
LSPAR:=./ls.mega.parameters

# accessibility graphs
%.acc.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) -fa $< > $@

# shrep graphs
%.shrep.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.shrep.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.shrep.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.shrep.gspan : %.fa %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -fasta $< > $@

# merge gspans
%.gspan : %.shrep.gspan %.acc.gspan
	$(MERGE_GSPAN) -shrep $*.shrep.gspan -acc $*.acc.gspan > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : bitsize=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

%.affy : %.shrep.gspan
	cat $< | grep '^t' | awk '{print $$5}' > $@

endif

.PHONY: all ls cv classstats clean distclean

# do predictions for all PROTEINS
all: $(PERF_FILES) results_aucpr.csv

# do parameter line search for all PROTEINS
ls : $(PARAM_FILES)

# do crossvalidation
cv : $(CV_FILES)

# generate staticstics on positive/negative composition
classstats : summary.cstats $(CSTAT_FILES)

# test various stuff
test: test_data_full_A.fa test_data_full_A.pred.fa \
	test_data_full_A.perf test_data_full_A.cstats test_data_full_A.param

test_data_full_A.fa :
	cp $(FA_DIR)/$@ $@

test_data_full_A.pred.fa :
	cp $(FA_DIR)/$@ $@

# keep fasta, predictions and results
clean:
	-rm -rf $(MODELS) log *.gspan *.threshold* *.model *.feature *.affy

# delete all files
distclean: clean
	-rm -rf *.param *.fa *.perf *.pred *.svrout *.ls.fa *.log results_aucpr.csv

%.cv : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv : %.feature %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

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
	$(LINESEARCH) -fa $< -param $(LSPAR) -mf Makefile -of $@ -bindir $(BINDIR) 2> >(tee $@.log >&2)
endif

# subset original fasta
%.ls.fa : %.fa
	cat $< | \
	$(FASTAPL) -e 'print ">", $$head, "\t", $$seq, "\n"' | \
	$(SHUF) -n $(LINESEARCH_INPUT_SIZE) | \
	$(PERL) -ane \
	'$$seq = pop @F; $$head = join(" ", @F); print $$head, "\n", $$seq, "\n";' > \
	$@

%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

%.svrout : %.model %.pred.feature
	time $(SVRPREDICT) $*.pred.feature $< $@

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
	$(PERF) -confusion < $@.threshold > $@
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
	$(PERL) -e 'print join("\t", "$(BASENAME)", "$(SET)", "$(LT)", "$(LN)", "$(HT)", "$(HN)"),"\n"' > $@

summary.cstats : $(CSTAT_FILES)
	( $(PERL) -e 'print join("\t", "protein", "set", "negative threshold", "negative instances", "positive threshold", "positive instances"),"\n"'; \
	cat $^ | sort -k1,2 ) > $@
