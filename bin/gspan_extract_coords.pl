#! /usr/bin/env perl

# extract offsets for calculating positionwise margins from gspan
# output format tab separated: graph type, sequence id, number of vertices with offset, offset

$seqid = -1;

while (<>) {
    if (/^t/) {
        $seqid++;
        /SEQLEN (\d+)/;
        $nvertices=$1;
        $offset=0;
        print join("\t", "seq", $seqid, $nvertices, $offset), "\n";
    }
    elsif (/^s/) {
        /WSTART (\d+)/;
        $wstart = $1;
        /WEND (\d+)/;
        $wend = $1;
        $offset = $wend - ($wstart - 1);
        /WSIZE (\d+)/;
        $nvertices=$1;
        print join("\t", "shrep", $seqid, $nvertices, $offset), "\n";
    }
}
