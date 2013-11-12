# sequence numbering starts at 0
BEGIN{
	seqid=-1;
}
# new graph collection -> new sequence
# vertex and nucleotide numbering starts at 0
/^t/{
	seqid++;
	vertex_id=0;
	nt_pos=0;
}
# new entry on viewpoint vertice (omitting relational vertices)
/^v/&&!/\^/&&!/#/{
	print seqid, vertex_id++, $3, nt_pos++;
}
# non-viewpoint  vertices only increment nucleotide position
/^V/&&!/\^/&&!/#/{
	nt_pos++
}
# restart sequence numbering on new subgraph
# (this does not work with shifted shreps!)
/^s/||/^u/||/^w/{
	nt_pos=0
}