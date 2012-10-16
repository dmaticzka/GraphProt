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
RBIN:=/usr/local/R/2.15.1-lx/bin/R --vanilla
SVRTRAIN:=/home/maticzkd/src/libsvm-3.12/svm-train -s 3 -t 0 -m $(SVR_CACHE)
SVRPREDICT:=/home/maticzkd/src/libsvm-3.12/svm-predict
PERF:=/home/maticzkd/src/stat/perf
LINESEARCH:=$(PERL) $(BINDIR)/lineSearch.pl
COMBINEFEATURES:=$(PERL) $(BINDIR)/combineFeatures.pl
SHUF:=/home/maticzkd/src/coreutils-8.15/src/shuf
FASTAPL:=$(PERL) /usr/local/user/RNAtools/fastapl
#FASTAPL:=/home/maticzkd/co/RNAtools/fastapl
FASTA2GSPAN:=$(PERL) /usr/local/user/RNAtools/fasta2shrep_gspan.pl
#FASTA2GSPAN:=/home/maticzkd/co/RNAtools/fasta2shrep_gspan.pl
SVMSGDNSPDK:=/home/maticzkd/local/svmsgdnspdk_120925/svmsgdnspdk
CREATE_EXTENDED_ACC_GRAPH:=$(PERL) $(BINDIR)/create_accgraph/createExtendedGraph.pl
MERGE_GSPAN:=$(PERL) $(BINDIR)/merge_gspan.pl
CAT_TABLES:=$(PERL) /home/maticzkd/co/MiscScripts/catTables.pl
FILTER_FEATURES:=$(PERL) $(BINDIR)/filter_features.pl
SUMMARIZE_MARGINS:=$(PERL) summarize_margins.pl
MARGINS2BG:=$(PERL) margins2bg.pl
BEDGRAPH2BIGWIG:=/usr/local/ucsctools/2012-02/bin/bedGraphToBigWig

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
# parameter files (from linesearch or default values)
PARAM_FILES:=$(patsubst %,%.train.param,$(BASENAMES))
# results of crossvalidations
CV_FILES:=$(patsubst %,%.train.cv,$(BASENAMES))
# models
MODEL_FILES:=$(patsubst %,%.train.model,$(BASENAMES))
# final results spearman correlation
CORRELATION_FILES:=$(patsubst %,%.test.correlation,$(BASENAMES))
# final results from perf
PERF_FILES:=$(patsubst %,%.test.perf,$(BASENAMES))
# nucleotide-wise margins
TESTPART_FILES:=$(patsubst %,%.test.nt_margins.summarized,$(BASENAMES))
# nucleotide-wise margins as bigWig
TESTPART_BIGWIG:=$(patsubst %,%.test.nt_margins.summarized.bw,$(BASENAMES))

## general feature and affinity creation (overridden where apropriate)
################################################################################
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy | %.param
	$(SVMSGDNSPDK) -a FEATURE -d $< -R $(RADIUS) -D $(DISTANCE) -gt $(DIRECTED)
	cat $<.feature | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $<.feature

# extract affinities from fasta
# expected to reside in last field of fasta header
%.affy : %.fa
	$(FASTAPL) -e '@ry = split(/\s/,$$head); print $$ry[-1], "\n"' < $< > $@
#	$(FASTAPL) -e 'print $$head[-1], "\n"' < $< > $@

## receipes specific to graph type
################################################################################
ifeq ($(GRAPH_TYPE),ONLYSEQ)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).onlyseq.parameters

%.gspan : %.fa
	$(FASTA2GSPAN) --seq-graph-t -nostr -stdout -fasta $< > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),STRUCTACC)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).structacc.parameters

%.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) -fa $< -W $(W_PRIMARY) -L $(L_PRIMARY) > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),SHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa | %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),PROBSHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa | %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) -stdout -q -Tp 0.05 -t $(ABSTRACTION) -M 5 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy | %.param
	# remove t and w, convert s to t
	cat $< | grep -v -e '^t' -e '^w' | sed 's/^s/t/' > $*_singleshreps
	# write out probabilities
	cat $< | grep '^s' | $(PERL) -ne '($$prob) = /SHAPEPROB ([0-9.]+)/; print $$prob, "\n"' > $*_probs
	# write out shrep membership
	cat $< | awk '/^t/{i++}/^s/{print i}' > $*_groups
	# compute features
	$(SVMSGDNSPDK) -a FEATURE -d $*_singleshreps -R $(RADIUS) -D $(DISTANCE) -gt $(DIRECTED)
	# compute probability-weighted features
	$(COMBINEFEATURES) $*_singleshreps.feature $*_probs $*_groups > $*
	# add affinities to features
	cat $* | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $*_singleshreps.feature
endif

################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).shrep_context.parameters

%.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan : %.fa | %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) -abstr -stdout -t $(ABSTRACTION) -M 5 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< > $@

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : bitsize=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : NHF=$(shell grep '^NHF ' $*.param | cut -f 2 -d' ')
%.feature : RR=$(shell grep '^RR ' $*.param | cut -f 2 -d' ')
%.feature : RD=$(shell grep '^RD ' $*.param | cut -f 2 -d' ')
%.feature : RW=$(shell grep '^RW ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan %.affy | %.param
	$(SVMSGDNSPDK) -kt ABSTRACT -a FEATURE -d $< -R $(RADIUS) -D $(DISTANCE) -gt $(DIRECTED) -anhf $(NHF) -rR $(RR) -rD $(RD) -rW $(RW)
	cat $<.feature | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $<.feature
endif

################################################################################
ifeq ($(GRAPH_TYPE),MEGA)
# line search parameters
LSPAR:=./ls.$(METHOD_ID).mega.parameters

# accessibility graphs
%.acc.gspan : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) -fa $< -W $(W_PRIMARY) -L $(L_PRIMARY) > $@

# shrep graphs
%.shrep.gspan : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.shrep.gspan : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.shrep.gspan : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.shrep.gspan : %.fa | %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) -stdout -t $(ABSTRACTION) -M 5 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< > $@

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
%.cv_svr : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.test.predictions_svr : %.train.model %.test.feature
	time $(SVRPREDICT) $*.test.feature $< $@

# affinities and predictions default format
%.predictions_affy : %.predictions_svr %.affy
	paste $*.affy $< > $@

# class membership and predictions default format
%.predictions_class : %.predictions_svr %.class
	paste $*.class $< > $@

endif


## support vector regression using sgd-derived subset of top features
################################################################################
ifeq ($(SVM),TOPSVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# train model; this one directly works on gspans
%.sgd_model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.sgd_model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.sgd_model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.sgd_model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.sgd_model : %.gspan %.class | %.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -a TRAIN -d $*.gspan -t $*.class -m $@ -R $(RADIUS) -D $(DISTANCE)

%.test.filter : %.train.filter
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
%.model : %.feature_filtered | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.test.predictions_svr : %.train.model %.test.feature_filtered
	time $(SVRPREDICT) $*.test.feature_filtered $< $@

# affinities and predictions default format
%.predictions_affy : %.predictions_svr %.affy
	# combine affinities and predictions
	paste $*.affy $< > $@

# class membership and predictions default format
%.predictions_class : %.predictions_svr %.class
	# combine affinities and predictions
	paste $*.class $< > $@
endif


## stochastic gradient descent
################################################################################
ifeq ($(SVM),SGD)
# extract single performance measure, used for linesearch decisions
%.cv : %.cv.perf
	cat $< | grep 'APR' | awk '{print $$NF}' > $@

# train model; this one directly works on gspans
%.model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.model : %.gspan %.class | %.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -a TRAIN -d $*.gspan -t $*.class -m $@ -R $(RADIUS) -D $(DISTANCE)

# evaluate model
%.test.predictions_sgd : RADIUS=$(shell grep '^R ' $*.test.param | cut -f 2 -d' ')
%.test.predictions_sgd : DISTANCE=$(shell grep '^D ' $*.test.param | cut -f 2 -d' ')
%.test.predictions_sgd : BITSIZE=$(shell grep '^b ' $*.test.param | cut -f 2 -d' ')
%.test.predictions_sgd : DIRECTED=$(shell grep '^DIRECTED ' $*.test.param | cut -f 2 -d' ')
%.test.predictions_sgd : %.train.model %.test.gspan %.test.class | %.test.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -R $(RADIUS) -D $(DISTANCE) -b $(BITSIZE) -a TEST -m $< -d $*.test.gspan -t $*.test.class
	mv $*.test.gspan.prediction $*.test.predictions_sgd

# affinities and predictions default format
%.predictions_affy : %.predictions_sgd %.affy
	cat $< | awk '{print $$2}' | paste $*.affy - > $@

# class membership and predictions default format: class{-1,1}, prediction
%.predictions_class : %.predictions_sgd %.class
	cat $< | awk '{print $$2}' | paste $*.class - > $@

# results from crossvalidation cast into default format: class{-1,1}, prediction
%.cv.predictions_class : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : %.gspan %.class | %.train.param
	time $(SVMSGDNSPDK) -gt $(DIRECTED) -b $(BITSIZE) -a CROSS_VALIDATION -cv $(CV_FOLD) -m $*.model -d $*.gspan -t $*.class -R $(RADIUS) -D $(DISTANCE) -sfx $*
	cat output.cv_predictions$* | awk '{print $$2==1?1:-1, $$4}' > $@
	-rm  -f output.cv_predictions$* $*.model_*

# compute margins of graph vertices
# vertex_margins format: seqid verticeid margin
%.test.vertex_margins : RADIUS=$(shell grep '^R ' $*.test.param | cut -f 2 -d' ')
%.test.vertex_margins : DISTANCE=$(shell grep '^D ' $*.test.param | cut -f 2 -d' ')
%.test.vertex_margins : BITSIZE=$(shell grep '^b ' $*.test.param | cut -f 2 -d' ')
%.test.vertex_margins : DIRECTED=$(shell grep '^DIRECTED ' $*.test.param | cut -f 2 -d' ')
%.test.vertex_margins : %.test.gspan %.test.class %.train.model | %.test.param
	$(SVMSGDNSPDK) -gt $(DIRECTED) -R $(RADIUS) -D $(DISTANCE) -b $(BITSIZE) -a TEST_PART -m $*.train.model -d $*.test.gspan -t $*.test.class
	mv $<.prediction_part $*.test.vertex_margins

# dictionary of all graph vertices
# dict file format: seqid v verticeid nt pos
%.vertex_dict : %.gspan
	cat $< | awk '/^t/{seqid++; vertex_id=0}/^v/{print seqid-1, vertex_id++, $$3, $$4}/^V/{print seqid-1, vertex_id++, $$3, $$4}' > $@

# compute nucleotide-wise margins from vertice margins
%.nt_margins : %.vertex_margins %.vertex_dict
	cat $< | $(PERL) ./vertex2ntmargins.pl -dict $*.vertex_dict | awk '$$2!=0' > $@

# format (tab separated): sequence id, sequence position, margin,
#                         min, max, mean, median, sum
%.nt_margins.summarized : %.nt_margins
	@echo ""
	@echo "summarizing nucleotide-wise margins:"
	$(SUMMARIZE_MARGINS) -W 25 < $< > $@

%.nt_margins.summarized.bg : %.nt_margins.summarized %.bed
	@echo ""
	@echo "converting margins to bedGraph"
	$(MARGINS2BG) -bed $*.bed < $< > $@

endif


## evaluations specific to RNAcompete analysis
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)

# class memberships {-1,0,1}
%.class : BASENAME=$(firstword $(subst _, ,$<))
%.class : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.class : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.class : %.affy
	cat $< | awk '{ if ($$1 > $(HT)) {print 1} else { if ($$1 < $(LT)) {print -1} else {print 0} } }' > $@

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
# combine input sequences
%.fa : %.positives.fa %.negatives.fa %.unknowns.fa
	( $(FASTAPL) -p -1 -e '$$head .= " 1";' < $*.positives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " -1";' < $*.negatives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " 0";' < $*.unknowns.fa ) > $@

%.bed : %.positives.bed %.negatives.bed %.unknowns.bed
	cat $^ > $@

%.positives.bed :
	@echo ""
	@echo "error: require file $@" && exit 255

%.unknowns.bed :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

%.negatives.bed :
	@echo ""
	@echo "using empty set of negatives!"
	touch $@

%.positives.fa :
	@echo ""
	@echo "error: require file $@" && exit 255

%.unknowns.fa :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

%.negatives.fa :
	@echo ""
	@echo "using empty set of negatives!"
	touch $@

# for clip data, affinities are actually the class
%.class : %.affy
	ln -sf $< $@
endif


## misc helper receipes
################################################################################
# we can save some disk space here
# %.gspan.gz : %.gspan
# 	gzip -f $<;

# # if available, create gspan from precomputed files
# %.gspan : %.gspan.gz
# 	zcat $< > $@

# link parameter file for simpler handling
%.test.param : %.train.param
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

# compute performance measures
# remove unknowns, set negative class to 0 for perf
%.perf : %.predictions_class
	cat $< | awk '$$1!=0' | sed 's/^-1/0/g' | $(PERF) -confusion > $@

# plot precision-recall
%.prplot : %.predictions_class
	cat $< | sed 's/^-1/0/g' | $(PERF) -plot pr | awk 'BEGIN{p=1}/ACC/{p=0}{if (p) {print}}' > $@

%.prplot.svg : %.prplot
	cat $< | gnuplot -e "set ylabel 'precision'; set xlabel 'recall'; set terminal svg; set style line 1 linecolor rgb 'black'; plot [0:1] [0:1] '-' using 1:2 with lines;" > $@

# compute correlation: correlation \t pvalue
%.correlation : %.predictions_affy
	cat $< | $(RBIN) --slave -e 'require(stats); data=read.table("$<", col.names=c("prediction","measurement")); t <- cor.test(data$$measurement, data$$prediction, method="spearman", alternative="greater"); write.table(cbind(t$$estimate, t$$p.value), file="$@", col.names=F, row.names=F, quote=F, sep="\t")'

results_aucpr.csv : $(PERF_FILES)
	grep -H -e APR -e ROC $^ | tr ':' "\t" | $(RBIN) --slave -e 'require(reshape); d<-read.table("stdin", col.names=c("id","variable","value")); write.table( cast(d), file="", row.names=F, quote=F, sep="\t")' > $@

results_correlation.csv : $(CORRELATION_FILES)
	$(CAT_TABLES) $(CORRELATION_FILES) > $@

# convert bedGraph to bigWig
%.bw : %.bg $(GENOME_BOUNDS)
	$(BEDGRAPH2BIGWIG) $*.bg $(GENOME_BOUNDS) $@

# do need genome bounds
$(GENOME_BOUNDS) :
	@echo ""
	@echo "error: require genome boundaries $@" && exit 255

## phony target section
################################################################################
.PHONY: all ls cv classstats test clean distclean

# do predictions and tests for all PROTEINS, summarize results
all: $(PERF_FILES) $(CORRELATION_FILES) results_aucpr.csv results_correlation.csv

# do parameter line search for all PROTEINS
ls : $(PARAM_FILES)

# do crossvalidation
cv : $(CV_FILES)

# train target
train : $(MODEL_FILES)

# test target
test : $(PERF_FILES) $(CORRELATION_FILES)

# compute nucleotide-wise margins
testpart : $(TESTPART_FILES)

# compute nucleotide-wise margins for genome-browser
testpart_coords : $(TESTPART_BIGWIG)

# generate staticstics on positive/negative composition
classstats : summary.cstats $(CSTAT_FILES)

# keep fasta, predictions and results
clean:
	-rm -rf log *.gspan *.threshold* *.feature *.affy *.feature_filtered \
	*.filter *.class

# delete all files
distclean: clean
	-rm -rf *.param *.perf *.predictions_class *.predictions_affy \
	*.predictions_svr *.predictions_sgd *.ls.fa *.log *.csv *model \
	*.sgeout *.class *.correlation *.cv *.cv.predictions \
	*.cv_svr *.model_* *.prplot *.prplot.svg

ifeq ($(EVAL_TYPE),CLIP)
# test various stuff
runtests: testclip.train.param testclip.train.cv \
	testclip.test.perf testclip.test.correlation \
	testclip.test.prplot.svg
endif
ifeq ($(EVAL_TYPE),RNACOMPETE)
# test various stuff
runtests: test_data_full_A.train.param test_data_full_A.train.cv \
	test_data_full_A.test.perf test_data_full_A.test.correlation \
	test_data_full_A.test.prplot.svg
endif

## insert additional rules into this file
################################################################################
include EXPERIMENT_SPECIFIC_RULES
