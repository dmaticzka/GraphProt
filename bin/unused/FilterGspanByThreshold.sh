#! /bin/bash 
[ $# -ne 4 ] &&  echo "Usage:  $(basename $0) <gspan filename> <Low threshold value> <High threshold value> <output basename>" && echo "Instead of: $(basename $0) $*"  && exit 1
datfn=$1
LT=$2
HT=$3
basename=$4

echo "FilterGspanByThreshold: processing $datfn"

# creating target file
zcat -f $datfn | grep "^t" | awk -v LT=$LT -v HT=$HT '{if($(NF)<LT){print $(NF)}else if($(NF)>HT){print $(NF)}else{print "nil"}}' > $basename.tmp.target
zcat -f $datfn | grep "^t" | awk -v LT=$LT -v HT=$HT '{if($(NF)<LT){print -1}else if($(NF)>HT){print 1}else{print 0}}' > $basename.tmp.class

# filtering out 0 targets
zcat -f $datfn | awk -v HF="$basename.tmp.target" 'BEGIN{ while (getline < HF ) {H[++c]=$1} } $1=="t"{n++} H[n]!="nil" {print $0}' | gzip > $basename.gspan.gz
zcat -f $basename.tmp.target | awk '$1!="nil"' > $basename.target
zcat -f $basename.tmp.class | awk '$1!=0' > $basename.class

# clean up
rm -rf $basename.tmp.target $basename.tmp.class
