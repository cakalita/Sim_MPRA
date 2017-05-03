#RNA (LCL)
dd <- unique(mpra_Filt[,c(1,13:18)])
DNA_sim <- simSimple1(prop=dd$DNA_prop,M=200,N=10000)
DNA_prop <- DNA_sim[,"R"]/(DNA_sim[,"R"]+DNA_sim[,"A"])
DNA_simm <- cbind(matrix(unlist(DNA_sim[,"R"]),ncol=1),matrix(unlist(DNA_sim[,"A"]),ncol=1))
DNA_sim_df <- data.frame(DNA_simm)
colnames(DNA_sim_df) <- c("DNA_R","DNA_A")
RNA_prop = DNA_prop
indNotNull <- sample(length(DNA_prop),round(0.001*length(DNA_prop)))
length(indNotNull)

deltaBeta <- rep(0,length(DNA_prop))
deltaBeta[indNotNull] = 2

RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
simRepAse201 <- simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/1+1))
simRepAse202 <- simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/2+1))
simRepAse203 <- simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/3+1))
simRepAse204 <- simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/4+1))
simRepAse205 <- simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/5+1))
simRepAse20df <- data.frame(cbind(simRepAse201,simRepAse202,simRepAse203,simRepAse204,simRepAse205))
colnames(simRepAse20df) <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
#dimnames(simRepAse20)[[3]] <- paste0("rep",1:5)
#simRepAse20m <- matrix(unlist(simRepAse20), ncol = 10, byrow = FALSE)
#simRepAse20df <- data.frame(simRepAse20m)
simRepAse20df$R_sum <- rowSums(simRepAse20df[,c(1,3,5,7,9)])
simRepAse20df$A_sum <- rowSums(simRepAse20df[,c(2,4,6,8,10)])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$R_sum,simRepAse20df_dna$A_sum,simRepAse20df_dna$DNA_prop)
#simttest_ase20 <- sim_ttest(x=simRepAse20df_dna)
#df2F <- subset(df2, R_sum>=5 & A_sum>=5)
table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
db <- cbind(simRepAse20df_dna, deltaBeta)
table(abs(db$deltaBeta)>0, fitQ$padj_quasar<0.1)

pval_ttest <- myttest(simRepAse20df_dna)
z <- cbind(simRepAse20df_dna,pval_ttest)
z$padj_ttest <- p.adjust(z$pval_ttest, method="BH")
all <- cbind(db, fitQ, z)

#pdf("qqplot_sim_LCL1_ASE20M10_db2_prop_0.001.pdf")
qq(all$pval3+1E-10,col ="white")
y=-log10((all$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=all$deltaBeta[o]*2+0.01)
y=-log10((all$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=all$deltaBeta[o]*2+0.01)
y <- sort(-log10( all$pval3[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='red' )
y <- sort(-log10( all$pval_ttest[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='purple')
legend("topleft", legend=c("QuASAR", "ttest","QuASAR null", "ttest null"),
    fill=c("green", "blue","red","purple"))
#dev.off()

#ind2
load("C:/Users/cakal/Google Drive/Lab/Tewhey/revisions/summed.RNA.LCL2_ttest.RData")
dd2 <- unique(mpra_2_Filt[,-c(2:7,10:11)])
DNA_sim <- simSimple1(prop=dd2$DNA_prop,M=200,N=10000)
DNA_prop <- DNA_sim[,"R"]/(DNA_sim[,"R"]+DNA_sim[,"A"])
DNA_simm <- cbind(matrix(unlist(DNA_sim[,"R"]),ncol=1),matrix(unlist(DNA_sim[,"A"]),ncol=1))
DNA_sim_df <- data.frame(DNA_simm)
colnames(DNA_sim_df) <- c("DNA_R","DNA_A")
RNA_prop = DNA_prop
indNotNull <- sample(length(DNA_prop),round(0.001*length(DNA_prop)))
length(indNotNull)

deltaBeta <- rep(0,length(DNA_prop))
deltaBeta[indNotNull] = 2

RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
simRepAse203 <- simSimple1(prop=RNA_prop,M=10,N=round((dd2$R+dd2$A)/3+1))
simRepAse204 <- simSimple1(prop=RNA_prop,M=10,N=round((dd2$R+dd2$A)/4+1))
simRepAse205 <- simSimple1(prop=RNA_prop,M=10,N=round((dd2$R+dd2$A)/5+1))
simRepAse20df <- data.frame(cbind(simRepAse203,simRepAse204,simRepAse205))
colnames(simRepAse20df) <- c("X1","X2","X3","X4","X5","X6")
#dimnames(simRepAse20)[[3]] <- paste0("rep",1:5)
#simRepAse20m <- matrix(unlist(simRepAse20), ncol = 10, byrow = FALSE)
#simRepAse20df <- data.frame(simRepAse20m)
simRepAse20df$R_sum <- rowSums(simRepAse20df[,c(1,3,5)])
simRepAse20df$A_sum <- rowSums(simRepAse20df[,c(2,4,6)])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$R_sum,simRepAse20df_dna$A_sum,simRepAse20df_dna$DNA_prop)
#simttest_ase20 <- sim_ttest(x=simRepAse20df_dna)
#df2F <- subset(df2, R_sum>=5 & A_sum>=5)
table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
db <- cbind(simRepAse20df_dna, deltaBeta)
table(abs(db$deltaBeta)>0, fitQ$padj_quasar<0.1)

myttest_ind2 <- function(simRep_null2_dna){
	rna_ref <- log2(simRep_null2_dna[,c(1,3,5)]/simRep_null2_dna[,"DNA_R"])
	rna_alt <- log2(simRep_null2_dna[,c(2,4,6)]/simRep_null2_dna[,"DNA_A"])
	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}
pval_ttest <- myttest(simRepAse20df_dna)
z <- cbind(simRepAse20df_dna,pval_ttest)
z$padj_ttest <- p.adjust(z$pval_ttest, method="BH")
all2 <- cbind(db, fitQ, z)

#pdf("qqplot_sim_LCL1_ASE20M10_db2_prop_0.001.pdf")
qq(all2$pval3+1E-10,col ="white")
y=-log10((all2$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=all2$deltaBeta[o]*2+0.01)
y=-log10((all2$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=all2$deltaBeta[o]*2+0.01)
y <- sort(-log10( all2$pval3[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='red' )
y <- sort(-log10( all2$pval_ttest[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='purple')
legend("topleft", legend=c("QuASAR", "ttest","QuASAR null", "ttest null"),
    fill=c("green", "blue","red","purple"))
#dev.off()

#combine
LCL1 <- cbind(all,all2)
combine <- function(x,y) return(pchisq(-2*log(x)-2*log(y),4,low=F))

LCL1 <- LCL[!is.na(LCL$betas_w.x), ]
LCL1 <- LCL1[!is.na(LCL1$betas_w.y), ]
LCL1 <- unique(LCL1)

#no DNA NA, removed
LCL <- transform(LCL1, DNA_prop=ifelse(DNA_prop.x=="NA", DNA_prop.y, DNA_prop.x), betas_T=((betas_w.x*betas.beta.binom.x) + (betas_w.y*betas.beta.binom.y))/(betas_w.x+betas_w.y), 
	betas_se_comb = sqrt(1/(betas_w.x+betas_w.y)))
LCL <- transform(LCL, betas_T_low = betas_T-(1.96*betas_se_comb), 
	betas_T_high = betas_T+(1.96*betas_se_comb),
	betas_z_comb = (betas_T - qlogis(DNA_prop)) /betas_se_comb)

LCL <- transform(LCL, p_comb = 2 * pnorm(-abs(betas_z_comb)))

LCL<- transform(LCL, comb_ttestpval=combine(pval_ttest,ttest_pval),comb_binompval=combine(pval_binom.x,pval_binom.y),
	comb_quasar_fisher=combine(pval3.x,pval3.y),comb_fisher_fisher=combine(pval_fisher.x,pval_fisher.y))

##############################################################
#############################################################3
myttest_real <- function(z){
  rna_ref <- log2(z[,c(3,5,7,9,11)]/z[,"DNA_R"])
  rna_alt <- log2(z[,c(2,4,6,8,10)]/z[,"DNA_A"])
  apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}
#z <- mpra_Filt[,c(1:2,13:14,17:18)]
mpra_Filt$ID <- paste0(mpra_Filt$rsID,"_",mpra_Filt$alt_hap)

lcl1 <- dcast(mpra_Filt, ID ~ Allele_class, value.var="NA12878_r1",fun.aggregate=sum)
lcl2 <- dcast(mpra_Filt, ID ~ Allele_class, value.var="NA12878_r2",fun.aggregate=sum)
lcl3 <- dcast(mpra_Filt, ID ~ Allele_class, value.var="NA12878_r3",fun.aggregate=sum)
lcl4 <- dcast(mpra_Filt, ID ~ Allele_class, value.var="NA12878_r4",fun.aggregate=sum)
lcl5 <- dcast(mpra_Filt, ID ~ Allele_class, value.var="NA12878_r5",fun.aggregate=sum)
counts <- cbind(lcl1,lcl2[,-1],lcl3[,-1],lcl4[,-1],lcl5[,-1])
counts1 <- merge(counts,unique(mpra_Filt[,c(13:15,22)]), by="ID")
pval_ttest <- myttest_real(counts1)
counts1$R_sum <- rowSums(counts1[,c(3,5,7,9,11)])
counts1$A_sum <- rowSums(counts1[,c(2,4,6,8,10)])
fitQ_real <- fitQuasarMpra(counts1$R_sum,counts1$A_sum,counts1$DNA_prop)

both <- cbind(counts1, fitQ_real, pval_ttest)
qq(both$pval3,col ="green")
y <- sort(-log10( mpra_Filt$pval_ttest))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')

load("C:/Users/cakal/Google Drive/Lab/Tewhey/revisions/summed.RNA.pvals.LCLcomb.RData")
qq(LCL$p_comb,col ="green")
y <- sort(-log10( LCL$comb_ttestpval))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')
