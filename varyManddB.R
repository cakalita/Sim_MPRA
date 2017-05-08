#RNA (LCL)
load("C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/summed.RNA.LCL1_ttest.RData")

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
simRepAse20 <- replicate(5,simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/5+1)))
simRepAse20df <- data.frame(simRepAse20)
simRepAse20df$r_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$a_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$r_sum,simRepAse20df_dna$a_sum,simRepAse20df_dna$DNA_prop)
table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
db <- cbind(simRepAse20df_dna, deltaBeta)

myttest <- function(simRep_null2_dna){
    rna_ref <- log2(simRep_null2_dna[ ,grepl("^R", names(simRep_null2_dna))]/simRep_null2_dna[,"DNA_R"])
    rna_alt <- log2(simRep_null2_dna[ ,grepl("^A", names(simRep_null2_dna))]/simRep_null2_dna[,"DNA_A"])
    apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}

pval_ttest_trial <- data.frame(pval_ttest=myttest_trial(simRepAse20df_dna))
pval_ttest_trial$padj_ttest <- p.adjust(pval_ttest_trial$pval_ttest, method="BH")

table(abs(db$deltaBeta)>0, pval_ttest_trial$padj_ttest<0.1)

M10_db2 <- cbind(db, fitQ, pval_ttest)
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
	file="C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/5_rep_varyManddB.RData")

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
simRepAse20df$r_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$a_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$r_sum,simRepAse20df_dna$a_sum,simRepAse20df_dna$DNA_prop)
db <- cbind(simRepAse20df_dna, deltaBeta)

#myttest_ind2 <- function(simRep_null2_dna){
#	rna_ref <- log2(simRep_null2_dna[,c(1,3,5)]/simRep_null2_dna[,"DNA_R"])
#	rna_alt <- log2(simRep_null2_dna[,c(2,4,6)]/simRep_null2_dna[,"DNA_A"])
#	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
#}
pval_ttest <- myttest(simRepAse20df_dna)
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

#RNA (LCL)
load("C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/summed.RNA.LCL1_ttest.RData")

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
deltaBeta[indNotNull] = 0.5
RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
simRepAse20 <- replicate(8,simSimple1(prop=RNA_prop,M=100,N=round((dd$R+dd$A)/5+1)))
simRepAse20df <- data.frame(simRepAse20)
simRepAse20df$r_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$a_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)
fitQ <- fitQuasarMpra(simRepAse20df_dna$r_sum,simRepAse20df_dna$a_sum,simRepAse20df_dna$DNA_prop)
db <- cbind(simRepAse20df_dna, deltaBeta)
pval_ttest_trial <- data.frame(pval_ttest=myttest_trial(simRepAse20df_dna))
pval_ttest_trial$padj_ttest <- p.adjust(pval_ttest_trial$pval_ttest, method="BH")

#M10_db2 <- cbind(db, fitQ, pval_ttest)
#M10_db1 <- cbind(db, fitQ, pval_ttest)
#M10_db0_5 <- cbind(db, fitQ, pval_ttest)
#M60_db2 <- cbind(db, fitQ, pval_ttest)
#M60_db1 <- cbind(db, fitQ, pval_ttest)
#M60_db0_5 <- cbind(db, fitQ, pval_ttest)
#M100_db2 <- cbind(db, fitQ, pval_ttest)
#M100_db1 <- cbind(db, fitQ, pval_ttest)
M100_db0_5 <- cbind(db, fitQ, pval_ttest)

save(M10_db2,M10_db1,M10_db0_5,M60_db2,M60_db1,M60_db0_5,
    M100_db2,M100_db1,M100_db0_5,
    file="C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/8_rep_varyManddB.RData")

M10_db2 <-  transform(M10_db2, M=10, dB=2)
M10_db1 <-  transform(M10_db1, M=10, dB=1)
M10_db0_5 <-  transform(M10_db0_5, M=10, dB=0.5)
M60_db2 <-  transform(M60_db2, M=60, dB=2)
M60_db1 <-  transform(M60_db1, M=60, dB=1)
M60_db0_5 <-  transform(M60_db0_5, M=60, dB=0.5)
M100_db2 <-  transform(M100_db2, M=100, dB=2)
M100_db1 <-  transform(M100_db1, M=100, dB=1)
M100_db0_5 <-  transform(M100_db0_5, M=100, dB=0.5)

#make list object of different simulations (changing variance and effect size (similar to fold change))
list8 <- list(M10_db2,M10_db1,M10_db0_5,M60_db2,M60_db1,M60_db0_5,M100_db2,M100_db1,M100_db0_5)

#set up series of padj thresholds
zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);

trial <-lapply(list8, function(i) {
    #i <- list1[[1]]
    NumTrue <- sum(abs(i$deltaBeta)>0, na.rm=T)
    aux <- sapply(zvec,function(z){
        NumDiscQ <- sum(i$padj_quasar<z,na.rm=T)
        NumDiscT <- sum(i$padj_ttest<z,na.rm=T) 
        NumTrueDiscQ <- sum(i$padj_quasar<z & abs(i$deltaBeta)>0, na.rm=T) 
        NumTrueDiscT <- sum(i$padj_ttest<z & abs(i$deltaBeta)>0, na.rm=T) 
        M <- i$M[1]
        deltaBeta = i$dB[1]
        data.frame(
            power_Q=NumTrueDiscQ/NumTrue,
            eFDR_Q=(NumDiscQ - NumTrueDiscQ)/NumDiscQ,
            power_T=NumTrueDiscT/NumTrue,
            eFDR_T=(NumDiscT - NumTrueDiscT)/NumDiscT,
            Z=z, 
            M=M,
            dB=deltaBeta) 
        })
    
})
zz8 <- lapply(trial,t)
zz8[] <- lapply(zz8,function(x) replace(x,is.na(x),0))

#plot ind 1 vs ind 2 dB1 M60
dat2 <- do.call("rbind", zz8)
zzm2 <- matrix(unlist(dat2),ncol=7,byrow = FALSE)
zzdf8 <- data.frame(zzm2)
colnames(zzdf8) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")

plot(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
    xlab="Z",ylab="power",main=paste0("M=60 dB=1 3,5,8 reps"))
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3)
points(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$power_Q[zzdf8$M==60 & zzdf8$dB==1], col='darkgreen',lty=2)
lines(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$power_Q[zzdf8$M==60 & zzdf8$dB==1], col='darkgreen',lty=2)
points(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$power_T[zzdf8$M==60 & zzdf8$dB==1], col='darkblue',lty=2)
lines(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$power_T[zzdf8$M==60 & zzdf8$dB==1], col='darkblue',lty=2)

plot(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
    xlab="Z",ylab="eFDR",main=paste0("M=60 dB=1 3,5,8 reps"))
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3)
points(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$eFDR_Q[zzdf8$M==60 & zzdf8$dB==1], col='darkgreen',lty=2)
lines(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$eFDR_Q[zzdf8$M==60 & zzdf8$dB==1], col='darkgreen',lty=2)
points(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$eFDR_T[zzdf8$M==60 & zzdf8$dB==1], col='darkblue',lty=2)
lines(x=zzdf8$Z[zzdf8$M==60 & zzdf8$dB==1], y=zzdf8$eFDR_T[zzdf8$M==60 & zzdf8$dB==1], col='darkblue',lty=2)
abline(0,1, col='red')

















#loop over set of simulations (possibly could work on this and add in the actual simulation with loop over M and dB conditions more automatically)
trial8_loops <-sapply(dBvec, function(d) {
    deltaBeta <- rep(0,length(DNA_prop))
    deltaBeta[indNotNull] = d
    RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
        aux <- sapply(Mvec, function(m, reps=8){
        simRepAse20 <- replicate(reps,simSimple1(prop=RNA_prop,M=m,N=round((dd$R+dd$A)/5+1)))
        simRepAse20df <- data.frame(simRepAse20)
        simRepAse20df$r_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
        simRepAse20df$a_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
        simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)
        fitQ <- fitQuasarMpra(simRepAse20df_dna$r_sum,simRepAse20df_dna$a_sum,simRepAse20df_dna$DNA_prop)
        #table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
        db <- cbind(simRepAse20df_dna, deltaBeta)
        pval_ttest_trial <- data.frame(pval_ttest=myttest_trial(simRepAse20df_dna))
        pval_ttest_trial$padj_ttest <- p.adjust(pval_ttest_trial$pval_ttest, method="BH")
        #i <- list8[[1]]
        NumTrue <- sum(abs(i$deltaBeta)>0, na.rm=T)
        M_dB <- cbind(db, fitQ, pval_ttest)
        M_dB <-  transform(M_dB, M=m, dB=d)
            aux2 <- sapply(zvec,function(z){
                NumDiscQ <- sum(M_dB$padj_quasar<z,na.rm=T)
                NumDiscT <- sum(M_dB$padj_ttest<z,na.rm=T) 
                NumTrueDiscQ <- sum(M_dB$padj_quasar<z & abs(M_dB$deltaBeta)>0, na.rm=T) 
                NumTrueDiscT <- sum(M_dB$padj_ttest<z & abs(M_dB$deltaBeta)>0, na.rm=T) 
                M <- M_dB$M[1]
                deltaBeta = M_dB$dB[1]
                data.frame(
                    power_Q=(NumTrueDiscQ/NumTrue),
                    eFDR_Q=(NumDiscQ - NumTrueDiscQ)/NumDiscQ,
                    power_T=(NumTrueDiscT/NumTrue),
                    eFDR_T=(NumDiscT - NumTrueDiscT)/NumDiscT,
                    Z=z, 
                    M=M,
                    dB=deltaBeta) 
            })
    })

})