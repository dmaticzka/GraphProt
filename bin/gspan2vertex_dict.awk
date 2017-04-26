# sequence numbering starts at 0
BEGIN{
	seqid=-1;
}
# new graph collection -> new sequence
# vertex and nucleotide numbering starts at 0
/^t/{
	seqid++;
	vertex_id=0;
}
# new entry on viewpoint vertice (omitting relational vertices)
/^v/ && !/\^/ && !/#/ {
	print seqid, vertex_id++, $3, $4-1;
}
# relational and abstract vertices increment vertex id and have no position
/^v/ && (/\^/ || /#/) {
	#print seqid, vertex_id++, $3, -1
	vertex_id++
}
# non-viewpoint vertices are silent
/^V/ && !/\^/ && !/#/ {
}
