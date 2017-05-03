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
simRepAse201 <- replicate(5,simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/5+1)))
simRepAse20df <- data.frame(simRepAse20)
simRepAse20df$R_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$A_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$R_sum,simRepAse20df_dna$A_sum,simRepAse20df_dna$DNA_prop)
table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
db <- cbind(simRepAse20df_dna, deltaBeta)
table(abs(db$deltaBeta)>0, fitQ$padj_quasar<0.1)

pval_ttest <- myttest(simRepAse20df_dna)
#pval_ttest$padj_ttest <- p.adjust(pval_ttest$pval_ttest, method="BH")
R3_M10_db2 <- cbind(db, fitQ, pval_ttest)
M10_db1 <- cbind(db, fitQ, pval_ttest)
M10_db0_5 <- cbind(db, fitQ, pval_ttest)
M10_db0_1 <- cbind(db, fitQ, pval_ttest)
M60_db2 <- cbind(db, fitQ, pval_ttest)
M60_db1 <- cbind(db, fitQ, pval_ttest)
M60_db0_5 <- cbind(db, fitQ, pval_ttest)
M60_db0_1 <- cbind(db, fitQ, pval_ttest)
M100_db2 <- cbind(db, fitQ, pval_ttest)
M100_db1 <- cbind(db, fitQ, pval_ttest)
M100_db0_5 <- cbind(db, fitQ, pval_ttest)
M100_db0_1 <- cbind(db, fitQ, pval_ttest)

save(M10_db2,M10_db1,M10_db0_5,M10_db0_1,M60_db2,M60_db1,M60_db0_5,M60_db0_1,
	M100_db2,M100_db1,M100_db0_5,M100_db0_1,
	file="5_rep_varyManddB.RData")

qq(M10_db2$pval_ttest+1E-10,col ="white")
y=-log10((M10_db2$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=M10_db2$deltaBeta[o]*2+0.01)
y=-log10((M10_db1$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=M10_db1$deltaBeta[o]*2+0.01)
y=-log10((M10_db0_5$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=M10_db0_5$deltaBeta[o]*2+0.01)
y=-log10((M10_db0_1$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=M10_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("M10_db2", "M10_db1","M10_db 0.5", "M10_db 0.1"),
    fill=c("green", "blue","red","purple"))

qq(M60_db2$pval_ttest+1E-10,col ="white")
y=-log10((M60_db2$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=M60_db2$deltaBeta[o]*2+0.01)
y=-log10((M60_db0_5$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=M60_db1$deltaBeta[o]*2+0.01)
y=-log10((M60_db0_5$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=M60_db0_5$deltaBeta[o]*2+0.01)
y=-log10((M60_db0_1$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=M60_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("M60_db2", "M60_db1","M60_db 0.5", "M60_db 0.1"),
    fill=c("green", "blue","red","purple"))

qq(M100_db2$pval_ttest+1E-10,col ="white")
y=-log10((M100_db2$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=M100_db2$deltaBeta[o]*2+0.01)
y=-log10((M100_db1$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=M100_db1$deltaBeta[o]*2+0.01)
y=-log10((M100_db0_5$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=M100_db0_5$deltaBeta[o]*2+0.01)
y=-log10((M100_db0_1$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=M100_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("M100_db2", "M100_db1","M100_db 0.5", "M100_db 0.1"),
    fill=c("green", "blue","red","purple"))


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
deltaBeta[indNotNull] = 0.1

RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
simRepAse20 <- replicate(3,simSimple1(prop=RNA_prop,M=10,N=round((dd2$R+dd2$A)/3+1)))
simRepAse20df <- data.frame(simRepAse20)
#simRepAse203 <- simSimple1(prop=RNA_prop,M=100,N=round((dd2$R+dd2$A)/3+1))
#simRepAse204 <- simSimple1(prop=RNA_prop,M=100,N=round((dd2$R+dd2$A)/3+1))
#simRepAse205 <- simSimple1(prop=RNA_prop,M=100,N=round((dd2$R+dd2$A)/3+1))
#simRepAse20df <- data.frame(cbind(simRepAse203,simRepAse204,simRepAse205))
simRepAse20df$R_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$A_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$R_sum,simRepAse20df_dna$A_sum,simRepAse20df_dna$DNA_prop)
db <- cbind(simRepAse20df_dna, deltaBeta)

myttest_ind2 <- function(simRep_null2_dna){
	rna_ref <- log2(simRep_null2_dna[,c(1,3,5)]/simRep_null2_dna[,"DNA_R"])
	rna_alt <- log2(simRep_null2_dna[,c(2,4,6)]/simRep_null2_dna[,"DNA_A"])
	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}
pval_ttest <- myttest_ind2(simRepAse20df_dna)
#pval_ttest$padj_ttest <- p.adjust(pval_ttest$pval_ttest, method="BH")

R3_M10_db2 <- cbind(db, fitQ, pval_ttest)
R3_M10_db1 <- cbind(db, fitQ, pval_ttest)
R3_M10_db0_5 <- cbind(db, fitQ, pval_ttest)
R3_M10_db0_1 <- cbind(db, fitQ, pval_ttest)
R3_M60_db2 <- cbind(db, fitQ, pval_ttest)
R3_M60_db1 <- cbind(db, fitQ, pval_ttest)
R3_M60_db0_5 <- cbind(db, fitQ, pval_ttest)
R3_M60_db0_1 <- cbind(db, fitQ, pval_ttest)
R3_M100_db2 <- cbind(db, fitQ, pval_ttest)
R3_M100_db1 <- cbind(db, fitQ, pval_ttest)
R3_M100_db0_5 <- cbind(db, fitQ, pval_ttest)
R3_M100_db0_1 <- cbind(db, fitQ, pval_ttest)

save(R3_M10_db2,R3_M10_db1,R3_M10_db0_5,R3_M10_db0_1,R3_M60_db2,R3_M60_db1,R3_M60_db0_5,R3_M60_db0_1,
	R3_M100_db2,R3_M100_db1,R3_M100_db0_5,R3_M100_db0_1,
	file="3_rep_varyManddB.RData")

qq(R3_M10_db2$pval3+1E-10,col ="white")
y=-log10((R3_M10_db2$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=R3_M10_db2$deltaBeta[o]*2+0.01)
y=-log10((R3_M10_db1$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=R3_M10_db1$deltaBeta[o]*2+0.01)
y=-log10((R3_M10_db0_5$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=R3_M10_db0_5$deltaBeta[o]*2+0.01)
y=-log10((R3_M10_db0_1$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=R3_M10_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("R3_M10_db2", "R3_M10_db1","R3_M10_db 0.5", "R3_M10_db 0.1"),
    fill=c("green", "blue","red","purple"))

qq(R3_M60_db2$pval3+1E-10,col ="white")
y=-log10((R3_M60_db2$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=R3_M60_db2$deltaBeta[o]*2+0.01)
y=-log10((R3_M60_db0_5$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=R3_M60_db1$deltaBeta[o]*2+0.01)
y=-log10((R3_M60_db0_5$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=R3_M60_db0_5$deltaBeta[o]*2+0.01)
y=-log10((R3_M60_db0_1$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=R3_M60_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("R3_M60_db2", "R3_M60_db1","R3_M60_db 0.5", "R3_M60_db 0.1"),
    fill=c("green", "blue","red","purple"))

qq(R3_M100_db2$pval3+1E-10,col ="white")
y=-log10((R3_M100_db2$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=R3_M100_db2$deltaBeta[o]*2+0.01)
y=-log10((R3_M100_db1$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=R3_M100_db1$deltaBeta[o]*2+0.01)
y=-log10((R3_M100_db0_5$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='red', cex=R3_M100_db0_5$deltaBeta[o]*2+0.01)
y=-log10((R3_M100_db0_1$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='purple', cex=R3_M100_db0_1$deltaBeta[o]*2+0.01)
legend("topleft", legend=c("R3_M100_db2", "R3_M100_db1","R3_M100_db 0.5", "R3_M100_db 0.1"),
    fill=c("green", "blue","red","purple"))

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

#combine going to work with M100_db0.5
M100_db1$betas_w <- 1/M100_db1$betas_se^2
R3_M100_db1$betas_w <- 1/R3_M100_db1$betas_se^2

LCL <- merge(M100_db1,R3_M100_db1, by="DNA_prop")
combine <- function(x,y) return(pchisq(-2*log(x)-2*log(y),4,low=F))

LCL1 <- LCL[!is.na(LCL$betas_w.x), ]
LCL1 <- LCL1[!is.na(LCL1$betas_w.y), ]
LCL1 <- unique(LCL1)

#no DNA NA, removed
LCL <- transform(LCL1, betas_T=((betas_w.x*betas.beta.binom.x) + (betas_w.y*betas.beta.binom.y))/(betas_w.x+betas_w.y), 
	betas_se_comb = sqrt(1/(betas_w.x+betas_w.y)))
LCL <- transform(LCL, betas_T_low = betas_T-(1.96*betas_se_comb), 
	betas_T_high = betas_T+(1.96*betas_se_comb),
	betas_z_comb = (betas_T - qlogis(DNA_prop)) /betas_se_comb)

LCL <- transform(LCL, p_comb = 2 * pnorm(-abs(betas_z_comb)),comb_ttestpval=combine(pval_ttest.x,pval_ttest.y))
#LCL<- transform(LCL, comb_ttestpval=combine(pval_ttest.x,pval_ttest.y),comb_binompval=combine(pval_binom.x,pval_binom.y),
	#comb_quasar_fisher=combine(pval3.x,pval3.y),comb_fisher_fisher=combine(pval_fisher.x,pval_fisher.y))

qq(LCL$p_comb,col ="green")
y <- sort(-log10( LCL$comb_ttestpval))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')
legend("topleft", legend=c("QuASAR", "ttest"), 
    fill=c("green", "blue"))