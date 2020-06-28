rm(list=ls())
load('./data/samp2pep.RData')
load('./data/prot-map.RData')
load('./data/peptideid.acc-df.RData')
load('./data/xx.ENSP.RData')
load('./data/xx.refseq.RData')
load('./data/xx.uniprot.df.RData')


map.pep2entrezID.study <- function(pep.obs,study.id){
	message("Processing ",study.id, " !")
	idx <- pepid.acc.df[,2] %in% pep.obs
	pep.acc<- pepid.acc.df[idx,1]
	# pepAcc to proteinAcc
	idx.pro <- prot_map[,1] %in% pep.acc
	prot.acc <- prot_map[idx.pro,]
	idx.ensp <- grep("^ENSP", prot.acc[,2])
	ensp <- unique(prot.acc[idx.ensp,2])
	ensp.pep <- unique(prot.acc[idx.ensp,1])
	getENSP.count <- function(x){
		# ENSP00000340878.4
		strsplit(x,split='\\.')[[1]][1]
	}
	ENSP <- sapply(ensp,getENSP.count)
	mapped.en <- xx.ENSP[names(xx.ENSP) %in% ENSP]
	entrezID.en <- unique(unlist(mapped.en))
	#----------------
	# REFSEQ
	idx.np <- grep("_", prot.acc[,2])
	np <- unique(prot.acc[idx.np,2])
	np.pep <- unique(prot.acc[idx.np,1])
	getNP.count <- function(x){
		# NP_001037.1
		strsplit(x,split='\\.')[[1]][1]
	}
	NP <- sapply(np,getNP.count)
	# REFSEQ TO ENTREZID
	mapped.np <- xx.refseq[names(xx.refseq) %in% NP]
	entrezID.np <- unique(unlist(mapped.np))
	#--------------------
	# UNIPROT
	up <- prot.acc[-c(idx.ensp,idx.np),2]
	up.pep <- unique(prot.acc[-c(idx.ensp,idx.np),1])
	getUP.count <- function(x){
		# NP_001037.1
		strsplit(x,split='\\-')[[1]][1]
	}
	UP <- unique(sapply(up,getUP.count))
	# UNIPROT TO ENTREZID
	mapped.up <- xx.uniprot.df[xx.uniprot.df$uniprot %in% UP,'entrezID']
	entrezID.up <- unique(mapped.up)
	#-----------------------------------
	# SUMMARY
	EN.c <- unique(entrezID.en)
	NP.c <- unique(entrezID.np)
	UP.c <- unique(entrezID.up)
	df.out <- data.frame(
		study.id = study.id,
		EN.peptide = length(ensp.pep),
		NP.peptide = length(np.pep),
		UP.peptide = length(up.pep),
		EN.protein = length(unique(ENSP)),
		NP.protein = length(unique(NP)),
		UP.protein = length(unique(UP)),
		EN.gene = length(EN.c),
		NP.gene = length(NP.c),
		UP.gene = length(UP.c),
		total.gene = length(unique(c(EN.c,NP.c,UP.c))),
		stringsAsFactors =FALSE
	)
	return(df.out)
}

#-----------------------------------
# SINGLE STUDYID
# pepID to pepAcc
library(doMC)
registerDoMC(10)
study.summary <- foreach(i = 1:length(samp2pep))%dopar%{
	pep.obs<-samp2pep[[i]]
	study.id <- names(samp2pep)[i]
	cat(i,'. ')
	map.pep2entrezID.study(pep.obs=pep.obs,study.id=study.id)
}
registerDoSEQ()
study.summary.df <- plyr::ldply(study.summary)
save(study.summary.df, file='./data/study.summary.df.RData')
