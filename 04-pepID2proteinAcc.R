rm(list=ls())
library(org.Hs.eg.db)
#-----------
# REFSEQ to EG
x.ENSP <- org.Hs.egENSEMBLPROT2EG
mapped_seqs <- mappedkeys(x.ENSP)
xx.ENSP <- as.list(x.ENSP[mapped_seqs])
save(xx.ENSP, file ='./data/xx.ENSP.RData')

#--------------------
# REFSEQ TO ENTREZID

x.refseq <- org.Hs.egREFSEQ2EG
mapped_seqs <- mappedkeys(x.refseq)
xx.refseq <- as.list(x.refseq[mapped_seqs])
save(xx.refseq, file ='./data/xx.refseq.RData')

#--------------------
# UNIPROT TO ENTREZID
x.uniprot <- org.Hs.egUNIPROT
mapped_seqs <- mappedkeys(x.uniprot)
xx.uniprot <- as.list(x.uniprot[mapped_seqs])
xx.len <- sapply(xx.uniprot,length)
xx.upChar <- rep(names(xx.uniprot),xx.len)
xx.up.df <- data.frame(
	uniprot = unlist(xx.uniprot),
	entrezID = xx.upChar,
	stringsAsFactors=FALSE
)
xx.uniprot.df <- xx.up.df
save(xx.uniprot.df, file ='./data/xx.uniprot.df.RData')

