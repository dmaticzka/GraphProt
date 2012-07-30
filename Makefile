## user defined  parameters are located here
################################################################################
include PARAMETERS


## general behaviour
################################################################################
SHELL:=/bin/bash
.DELETE_ON_ERROR:
ifeq ($(SECONDARY),YES)
# don't delete intermediate files
.SECONDARY:
endif


## paths
################################################################################
PROJDIR:=/home/maticzkd/projects/RBPaffinity
FA_DIR:=$(PROJDIR)/data/fasta
THR_DIR:=$(PROJDIR)/data/thresholds/
# expect binaries to reside in pwd, otherwise this variable must be overwritten
BINDIR:=$(shell pwd)


## binaries
################################################################################
PERL:=/usr/local/perl/bin/perl
RBIN:=/usr/local/R/2.12.1-lx/bin/R
SVRTRAIN:=/home/maticzkd/src/libsvm-3.12/svm-train -s 3 -t 0 -m $(SVR_CACHE)
SVRPREDICT:=/home/maticzkd/src/libsvm-3.12/svm-predict
PERF:=/home/maticzkd/src/stat/perf
LINESEARCH:=$(PERL) $(BINDIR)/lineSearch.pl
COMBINEFEATURES:=$(PERL) $(BINDIR)/combineFeatures.pl
SHUF:=/home/maticzkd/src/coreutils-8.15/src/shuf
FASTAPL:=$(PERL) /usr/local/user/RNAtools/fastapl
#FASTAPL:=/home/maticzkd/repositories/RNAtools/fastapl
FASTA2GSPAN:=$(PERL) /usr/local/user/RNAtools/fasta2shrep_gspan.pl
#FASTA2GSPAN:=/home/maticzkd/repositories/RNAtools/fasta2shrep_gspan.pl
SVMSGDNSPDK:=/home/maticzkd/repositories/svmsgdnspdk_dir_nspdkvp/svmsgdnspdk
CREATE_EXTENDED_ACC_GRAPH:=$(PERL) $(BINDIR)/create_accgraph/createExtendedGraph.pl
MERGE_GSPAN:=$(PERL) $(BINDIR)/merge_gspan.pl
CAT_TABLES:=$(PERL) /home/maticzkd/repositories/MiscScripts/catTables.pl
FILTER_FEATURES:=$(PERL) $(BINDIR)/filter_features.pl


## set appropriate id (used to determine which parameter sets to use)
################################################################################
ifeq ($(SVM),SVR)
METHOD_ID=svr
endif
ifeq ($(SVM),TOPSVR)
METHOD_ID=svr
endif
ifeq ($(SVM),SGD)
METHOD_ID=sgd
endif


## set targets for RNAcompete evaluation
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)
# filenames for full data sets
FULL_BASENAMES:=$(patsubst %,%_data_full_A,$(PROTEINS)) \
			$(patsubst %,%_data_full_B,$(PROTEINS))

# filenames of data sets containing only weakly structured sequences
BRUIJN_BASENAMES:=$(patsubst %,%_data_bruijn_A,$(PROTEINS)) \
			$(patsubst %,%_data_bruijn_B,$(PROTEINS))

# extract prefixes for further assembling of target filenames
ifeq ($(TRAINING_SETS),FULL)
BASENAMES:=$(FULL_BASENAMES)
else
ifeq ($(TRAINING_SETS),WEAK)
BASENAMES:=$(BRUIJN_BASENAMES)
else
BASENAMES:=$(FULL_BASENAMES) $(BRUIJN_BASENAMES)
endif
endif

# general class statistics
CSTAT_FILES:=$(patsubst %,%.cstats,$(FULL_BASENAMES))
endif


## set targets for CLIP-seq evaluation
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
BASENAMES=$(PROTEINS)
endif


## set targets common to all evaluations
################################################################################
# final results spearman correlation
CORRELATION_FILES:=$(patsubst %,%.correlation,$(BASENAMES))
# final results from perf
PERF_FILES:=$(patsubst %,%.perf,$(BASENAMES))
# parameter files (from linesearch or default values)
PARAM_FILES:=$(patsubst %,%.param,$(BASENAMES))
# results of crossvalidations
CV_FILES:=$(patsubst %,%.cv,$(BASENAMES))


## general feature and affinity creation (overridden where apropriate)
################################################################################
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy %.param
	$(SVMSGDNSPDK) -a FEATUREGENERATION -d $< -ll 1 $(RADIUS) $(DISTANCE) -gt $(DIRECTED) -pfx $*
	cat R$(RADIUS)D$(DISTANCE)$*output.vec | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f R$(RADIUS)D$(DISTANCE)$*output.vec

# extract affinities from fasta
# expected to reside in last field of fasta header
%.affy : %.fa
	$(FASTAPL) -e 'print $$head[-1], "\n"' < $< > $@

# combine input sequences
%.fa : %.positives.fa %.negatives.fa %.unknowns.fa
	( $(FASTAPL) -p -1 -e '$$head .= " 1";' < $*.positives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " -1";' < $*.negatives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " 0";' < $*.unknowns.fa ) > $@

%.unknowns.fa :
	echo "doing supervised training only"
	touch $@

## receipes specific to graph type
################################################################################
ifeq ($(GRAPH_TYPE),ONLYSEQ)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).onlyseq.parameters

%.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) --nostruct -fa $< > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),STRUCTACC)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).structacc.parameters

%.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) -fa $< > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),SHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -fasta $< > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),PROBSHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep.parameters

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
endif

################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep_context.parameters

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
endif

################################################################################
ifeq ($(GRAPH_TYPE),MEGA)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).mega.parameters

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
endif


## receipes specific to SVM type
################################################################################
# support vector regression
ifeq ($(SVM),SVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.svrout : %.model %.pred.feature
	time $(SVRPREDICT) $*.pred.feature $< $@

# affinities and predictions: default format
%.pred : %.svrout %.pred.affy
	# combine affinities and predictions
	paste $*.pred.affy $< > $@
endif


## stochastic gradient descent
################################################################################
ifeq ($(SVM),SGD)
# results from crossvalidation
%.cv_sgd : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_sgd : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_sgd : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.cv_sgd : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.cv_sgd : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.cv_sgd : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.cv_sgd : %.gspan %.class %.param
	time $(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -mode FILE -a CROSSVALIDATION -cv $(CV_FOLD) -d $*.gspan -t $*.class -ll 1 $(RADIUS) $(DISTANCE) -pfx $*
	cat $*output.cv_predictions |awk '{print $$2==1?1:0, $$4}' | $(PERF) -confusion > $@
	-rm -f $*output.cv_predictions $*model_*

%.cv : %.cv_sgd
	cat $< | grep 'APR' | awk '{print $$NF}' > $@

# train model; this one directly works on gspans
%.model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.model : %.gspan %.class %.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -mode FILE -a TRAIN -d $*.gspan -t $*.class -m $@ -ll 1 $(RADIUS) $(DISTANCE)

# this version of SGD reads all parameters from model
%.output.predictions : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.output.predictions : %.model %.pred.gspan %.pred.class
	$(SVMSGDNSPDK) -gt $(DIRECTED) -mode FILE -a TEST -m $< -d $*.pred.gspan -t $*.pred.class -pfx $*.

# affinities and predictions: default format
%.pred : %.output.predictions %.pred.affy
	cat $< | awk '{print $$2}' | paste $*.pred.affy - > $@
endif


## support vector regression using sgd-derived subset of top features
################################################################################
ifeq ($(SVM),TOPSVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# train model; this one directly works on gspans
%.sgd_model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.sgd_model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.sgd_model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.sgd_model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.sgd_model : %.gspan %.class %.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -mode FILE -a TRAIN -d $*.gspan -t $*.class -m $@ -ll 1 $(RADIUS) $(DISTANCE)

%.pred.filter : %.filter
	ln -s $< $@

%.filter : NFEAT=$(shell cat $< | grep '^w ' | sed 's/^w //' | tr ' :' "\n\t" | wc -l)
%.filter : TENP=$(shell echo "$(NFEAT) / 5" | bc)
%.filter : %.sgd_model
	cat $< | grep '^w ' | sed 's/^w //' | tr ' :' "\n\t" | sort -k2,2gr | head -n $(TENP) | cut -f 1 | sort -n > $@

%.feature_filtered : %.feature %.filter
	$(FILTER_FEATURES) --feature $< --filter $*.filter > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature_filtered %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.svrout : %.model %.pred.feature_filtered
	time $(SVRPREDICT) $*.pred.feature_filtered $< $@

# affinities and predictions: default format
%.pred : %.svrout %.pred.affy
	# combine affinities and predictions
	paste $*.pred.affy $< > $@
endif


## evaluations specific to RNAcompete analysis
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)
# test various stuff
test: test_data_full_A.fa test_data_full_A.pred.fa \
	test_data_full_A.perf test_data_full_A.correlation \
	test_data_full_A.cstats test_data_full_A.param

# helper receipes for test
test_data_full_A.fa :
	cp -f $(FA_DIR)/$@ $@
test_data_full_A.pred.fa :
	cp -f $(FA_DIR)/$@ $@
test_data_full_B.fa :
	cp -f $(FA_DIR)/$@ $@
test_data_full_B.pred.fa :
	cp -f $(FA_DIR)/$@ $@

# class memberships {-1,0,1}
%.class : BASENAME=$(firstword $(subst _, ,$<))
%.class : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.class : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.class : %.affy
	cat $< | awk '{ if ($$1 > $(HT)) {print 1} else { if ($$1 < $(LT)) {print -1} else {print 0} } }' > $@

# compute performance measures
# .perf is tab separated: affinity \t prediction
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

# some statistics about class distribution
%.cstats : BASENAME=$(firstword $(subst _, ,$<))
%.cstats : TYPE=$(word 3,$(subst _, ,$<))
%.cstats : SET=$(word 4,$(subst ., ,$(subst _, ,$<)))
%.cstats : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.cstats : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.cstats : HN=$(shell cat $< | grep '^>' | awk '$$NF > $(HT)' | wc -l)
%.cstats : LN=$(shell cat $< | grep '^>' | awk '$$NF < $(LT)' | wc -l)
%.cstats : %.fa
	$(PERL) -e 'print join("\t", "$(BASENAME)", "$(SET)", "$(LT)", "$(LN)", "$(HT)", "$(HN)"),"\n"' > $@

# final class summary
summary.cstats : $(CSTAT_FILES)
	( $(PERL) -e 'print join("\t", "protein", "set", "negative threshold", "negative instances", "positive threshold", "positive instances"),"\n"'; \
	cat $^ | sort -k1,2 ) > $@
endif


## evaluations specific to CLIP analysis
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
# test various stuff
test: testclip.fa testclip.gspan testclip.affy testclip.class testclip.cv

# for clip data, affinities are actually the class
%.class : %.affy
	ln -sf $< $@
endif


## misc helper receipes
################################################################################
# we can save some disk space here
# %.gspan.gz : %.gspan
# 	gzip -f $<;

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

# subset fastas prior to line search
%.ls.fa : %.fa
	cat $< | \
	$(FASTAPL) -e 'print ">", $$head, "\t", $$seq, "\n"' | \
	$(SHUF) -n $(LINESEARCH_INPUT_SIZE) | \
	$(PERL) -ane \
	'$$seq = pop @F; $$head = join(" ", @F); print $$head, "\n", $$seq, "\n";' > \
	$@

# final results summary
%.correlation : %.pred
	cat $< | $(RBIN) --slave -e 'data=read.table("$<", col.names=c("prediction","measurement")); t <- cor.test(data$$measurement, data$$prediction, method="spearman", alternative="greater"); write.table(cbind(t$$estimate, t$$p.value), file="$@", col.names=F, row.names=F, quote=F, sep="\t")'

results_aucpr.csv : $(PERF_FILES)
	grep -H -e APR -e ROC $^ | tr ':' "\t" | $(RBIN) --slave -e 'require(reshape); d<-read.table("stdin", col.names=c("id","variable","value")); write.table( cast(d), file="", row.names=F, quote=F, sep="\t")' > $@

results_correlation.csv : $(CORRELATION_FILES)
	$(CAT_TABLES) $(CORRELATION_FILES) > $@


## phony target section
################################################################################
.PHONY: all ls cv classstats test clean distclean

# do predictions for all PROTEINS
all: $(PERF_FILES) $(CORRELATION_FILES) results_aucpr.csv results_correlation.csv

# do parameter line search for all PROTEINS
ls : $(PARAM_FILES)

# do crossvalidation
cv : $(CV_FILES)

# generate staticstics on positive/negative composition
classstats : summary.cstats $(CSTAT_FILES)

# keep fasta, predictions and results
clean:
	-rm -rf log *.gspan *.threshold* *.feature *.affy *.feature_filtered \
	*.filter

# delete all files
distclean: clean
	-rm -rf *.param *.fa *.perf *.pred *.svrout *.ls.fa *.log *.csv *model \
	*.sgeout *.class *.output.predictions *.correlation *.cv *.cv_sgd \
	*.cv_svr
