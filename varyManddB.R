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

#pval_ttest_trial <- data.frame(pval_ttest=myttest_trial(simRepAse20df_dna))
#pval_ttest_trial$padj_ttest <- p.adjust(pval_ttest_trial$pval_ttest, method="BH")

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
deltaBeta[indNotNull] = 1
RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))
simRepAse20 <- replicate(8,simSimple1(prop=RNA_prop,M=60,N=round((dd$R+dd$A)/5+1)))
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
#M60_db2 <- cbind(db, fitQ, pval_ttest_trial)
M60_db1 <- cbind(db, fitQ, pval_ttest_trial)
#M60_db0_5 <- cbind(db, fitQ, pval_ttest_trial)
#M100_db2 <- cbind(db, fitQ, pval_ttest)
#M100_db1 <- cbind(db, fitQ, pval_ttest)
#M100_db0_5 <- cbind(db, fitQ, pval_ttest)

#save(M10_db2,M10_db1,M10_db0_5,M60_db2,M60_db1,M60_db0_5,
#    M100_db2,M100_db1,M100_db0_5,
#    file="C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/8_rep_varyManddB.RData")
save(M60_db2,M60_db1,M60_db0_5,
    file="C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/8_rep_2nd_varyManddB.RData")


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
list8_2 <- list(M60_db2,M60_db1,M60_db0_5)

#set up series of padj thresholds
zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);

trial <-lapply(list8_2, function(i) {
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
rowCdf <- data.frame(zzm2)
colnames(rowCdf) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")

plot(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$M==60 & rowCdf$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
    xlab="Z",ylab="power",main=paste0("M=60 dB=1 3,5,8 reps"))
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$M==60 & rowCdf$dB==1], col='lightgreen',lty=3)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$M==60 & rowCdf$dB==1], col='lightblue',lty=3)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$M==60 & rowCdf$dB==1], col='lightblue',lty=3)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green',lty=1)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green',lty=1)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='blue',lty=1)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='blue',lty=1)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$M==60 & rowCdf$dB==1], col='darkgreen',lty=2)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$M==60 & rowCdf$dB==1], col='darkgreen',lty=2)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$M==60 & rowCdf$dB==1], col='darkblue',lty=2)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$M==60 & rowCdf$dB==1], col='darkblue',lty=2)

plot(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$M==60 & rowCdf$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
    xlab="Z",ylab="eFDR",main=paste0("M=60 dB=1 3,5,8 reps"))
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$M==60 & rowCdf$dB==1], col='lightgreen',lty=3)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$M==60 & rowCdf$dB==1], col='lightblue',lty=3)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$M==60 & rowCdf$dB==1], col='lightblue',lty=3)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green',lty=1)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green',lty=1)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='blue',lty=1)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='blue',lty=1)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$M==60 & rowCdf$dB==1], col='darkgreen',lty=2)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$M==60 & rowCdf$dB==1], col='darkgreen',lty=2)
points(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$M==60 & rowCdf$dB==1], col='darkblue',lty=2)
lines(x=rowCdf$Z[rowCdf$M==60 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$M==60 & rowCdf$dB==1], col='darkblue',lty=2)
abline(0,1, col='red')















dBvec <- c(0.5,1,2)
Mvec <- c(10,60,100)
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
        M_dB <- cbind(db, fitQ, pval_ttest_trial)
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
            aux2 <- do.call(rbind,aux2)
            aux2
    })
    aux <- do.call(rbind,aux)
    aux
})
trial8_loops <- do.call(rbind,trial8_loops)


simulAndEval <- function(zvec,DNA_prop,m=60,d=1,reps=5,reads.div=5){
   deltaBeta <- rep(0,length(DNA_prop))
    deltaBeta[indNotNull] = d
    RNA_prop = plogis((qlogis(DNA_prop)+deltaBeta))

        simRepAse20 <- replicate(reps,simSimple1(prop=RNA_prop,M=m,N=round((dd$R+dd$A)/reads.div+1)))
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
        NumTrue <- sum(abs(deltaBeta)>0, na.rm=T)
        M_dB <- cbind(db, fitQ, pval_ttest_trial)
        M_dB <-  transform(M_dB, M=m, dB=d)
            aux2 <- sapply(zvec,function(z){
                NumDiscQ <- sum(M_dB$padj_quasar<z,na.rm=T)
                NumDiscT <- sum(M_dB$padj_ttest<z,na.rm=T) 
                NumTrueDiscQ <- sum(M_dB$padj_quasar<z & abs(M_dB$deltaBeta)>0, na.rm=T) 
                NumTrueDiscT <- sum(M_dB$padj_ttest<z & abs(M_dB$deltaBeta)>0, na.rm=T) 
                M <- M_dB$M[1]
                deltaBeta = M_dB$dB[1]
                c(
                    power_Q=(NumTrueDiscQ/NumTrue),
                    eFDR_Q=(NumDiscQ - NumTrueDiscQ)/NumDiscQ,
                    power_T=(NumTrueDiscT/NumTrue),
                    eFDR_T=(NumDiscT - NumTrueDiscT)/NumDiscT,
                    Z=z, 
                    M=M,
                    dB=deltaBeta,
                    reps=reps,
                    reads.div=reads.div
                    ) 
            })
            aux2 <- as.matrix(t(aux2))
            aux2 <- replace(aux2,is.na(aux2),0)
            aux2
}

aux2 <- simulAndEval(zvec,DNA_prop)


simulAndEvalRep <- function(zvec,DNA_prop,m=60,d=1,reps=5,reads.div=5,simReps=5){
    aux<- replicate(simReps,simulAndEval(zvec,DNA_prop,m=m,d=d,reps=reps,reads.div=reads.div))
    apply(aux,c(1,2),mean)
    #lapply(auxlist,function(x) replace(x,is.na(x),0))
    ##do.call(mean,auxlist)
    ##auxlist
}
#zz_div2[] <- lapply(zz_div2,function(x) replace(x,is.na(x),0))

aux2 <- simulAndEvalRep(zvec,DNA_prop,m=60,d=1,reps=5,reads.div=5,simReps=5)

rowA <- lapply(list(10,60,100),function(m){
    simulAndEvalRep(zvec,DNA_prop,m=m,d=1,reps=5,reads.div=5,simReps=5)
})
dat <- do.call("rbind", rowA)
rowAdf <- data.frame(dat)

rowB <- lapply(list(0.5,2),function(d){
    simulAndEvalRep(zvec,DNA_prop,m=60,d=d,reps=5,reads.div=5,simReps=5)
})
dat <- do.call("rbind", rowB)
rowBdf <- data.frame(dat)

rowC <- lapply(list(3,8),function(reps){
    simulAndEvalRep(zvec,DNA_prop,m=60,d=1,reps=reps,reads.div=5,simReps=5)
})
dat <- do.call("rbind", rowC)
rowCdf <- data.frame(dat)

rowD <- lapply(list(2,10),function(reads.div){
    simulAndEvalRep(zvec,DNA_prop,m=60,d=1,reps=5,reads.div=reads.div,simReps=5)
})
dat <- do.call("rbind", rowD)
rowDdf <- data.frame(dat)

save(rowAdf, rowBdf, rowCdf, rowDdf, file="simulations_5reps.RData")

png(paste0("Mvaried_dB1_zxPower_combineddB_noM100_linethick_nomargin.png"),width=800*4,height=900*4,pointsize=12*4*1.5)
#par(mfrow=c(5,3))
par(mar = (c(0, 0, 0, 0)+0.25), oma = c(4, 1, 3, 3))
layout(matrix(1:12,nrow=4, ncol=3, byrow=TRUE), widths=c(1.25,2,2))
    legText <- c("QuASAR-MPRA","M=10","M=60","M=100","T-test","M=10","M=60","M=100")
    # create a blank graph -- automatically scales -1 to +1 on both axes 
    barplot(0,0, axes=FALSE) 
    legend(x=0.8, y=0,
        title=paste0("M=varied dB=1"), 
           legend=legText, 
           col=c(NA,"green1","green3","green4",NA,"royalblue1","royalblue3","royalblue4"),
           lwd=c(NA,5*4,3*4,2*4,NA,5*4,3*4,2*4),
           xjust=1, 
           yjust=0.5) 
plot(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==10 & rowAdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1),
  xlab="",ylab="",xaxt='n',yaxt='n')
mtext("Power",side=3,line=1)
lines(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==10 & rowAdf$dB==1], col='green1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==10 & rowAdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==10 & rowAdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==100 & rowAdf$dB==1], col='green4',lwd=2*4)
lines(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==100 & rowAdf$dB==1], col='green4',lwd=2*4)
points(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==100 & rowAdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==100 & rowAdf$dB==1], col='royalblue4',lwd=2*4)
plot(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==10 & rowAdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1),
  xlab="",ylab="",xaxt='n',yaxt='n')
mtext("eFDR",side=3,line=1)
axis(4, las=1)
lines(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==10 & rowAdf$dB==1], col='green1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==10 & rowAdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowAdf$Z[rowAdf$M==10 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==10 & rowAdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==100 & rowAdf$dB==1], col='green4',lwd=2*4)
lines(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==100 & rowAdf$dB==1], col='green4',lwd=2*4)
points(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==100 & rowAdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowAdf$Z[rowAdf$M==100 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==100 & rowAdf$dB==1], col='royalblue4',lwd=2*4)
abline(0,1, col='red')
    legText <- c("QuASAR-MPRA","dB=0.5","dB=1","dB=2","T-test","dB=0.5","dB=1","dB=2")
    # create a blank graph -- automatically scales -1 to +1 on both axes 
    barplot(0,0, axes=FALSE) 
    legend(x=0.8, y=0, 
        title=paste0("M=60 dB=varied"),
           legend=legText, 
           col=c(NA,"green1","green3","green4",NA,"royalblue1","royalblue3","royalblue4"),
           lwd=c(NA,5*4,3*4,2*4,NA,5*4,3*4,2*4),
           xjust=1, 
           yjust=0.5) 
plot(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$power_Q[rowBdf$M==60 & rowBdf$dB==0.5], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1)
  ,xlab="",ylab="",xaxt='n',yaxt='n')
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$power_Q[rowBdf$M==60 & rowBdf$dB==0.5], col='green1',lwd=5*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$power_T[rowBdf$M==60 & rowBdf$dB==0.5], col='royalblue1',lwd=5*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$power_T[rowBdf$M==60 & rowBdf$dB==0.5], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$power_Q[rowBdf$M==60 & rowBdf$dB==2], col='green4',lwd=2*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$power_Q[rowBdf$M==60 & rowBdf$dB==2], col='green4',lwd=2*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$power_T[rowBdf$M==60 & rowBdf$dB==2], col='royalblue4',lwd=2*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$power_T[rowBdf$M==60 & rowBdf$dB==2], col='royalblue4',lwd=2*4)
plot(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$eFDR_Q[rowBdf$M==60 & rowBdf$dB==0.5], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1)
  ,xlab="",ylab="",xaxt='n',yaxt='n')
axis(4, las=1)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$eFDR_Q[rowBdf$M==60 & rowBdf$dB==0.5], col='green1',lwd=5*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$eFDR_T[rowBdf$M==60 & rowBdf$dB==0.5], col='royalblue1',lwd=5*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==0.5], y=rowBdf$eFDR_T[rowBdf$M==60 & rowBdf$dB==0.5], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$eFDR_Q[rowBdf$M==60 & rowBdf$dB==2], col='green4',lwd=2*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$eFDR_Q[rowBdf$M==60 & rowBdf$dB==2], col='green4',lwd=2*4)
points(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$eFDR_T[rowBdf$M==60 & rowBdf$dB==2], col='royalblue4',lwd=2*4)
lines(x=rowBdf$Z[rowBdf$M==60 & rowBdf$dB==2], y=rowBdf$eFDR_T[rowBdf$M==60 & rowBdf$dB==2], col='royalblue4',lwd=2*4)
abline(0,1, col='red')
    legText <- c("QuASAR-MPRA","3 reps","5 reps","8 reps","T-test","3 reps","5 reps","8 reps")
    # create a blank graph -- automatically scales -1 to +1 on both axes 
    barplot(0,0, axes=FALSE) 
    legend(x=0.8, y=0, 
      title=paste0("M=60 dB=1"),
           legend=legText, 
           col=c(NA,"green1","green3","green4",NA,"royalblue1","royalblue3","royalblue4"),
           lwd=c(NA,5*4,3*4,2*4,NA,5*4,3*4,2*4),
           xjust=1, 
           yjust=0.5) 
plot(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$reps==3 & rowCdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1)
  ,xlab="",ylab="",xaxt='n',yaxt='n')
lines(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$reps==3 & rowCdf$dB==1], col='green1',lwd=5*4)
points(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$reps==3 & rowCdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$reps==3 & rowCdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$reps==8 & rowCdf$dB==1], col='green4',lwd=2*4)
lines(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$power_Q[rowCdf$reps==8 & rowCdf$dB==1], col='green4',lwd=2*4)
points(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$reps==8 & rowCdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$power_T[rowCdf$reps==8 & rowCdf$dB==1], col='royalblue4',lwd=2*4)
plot(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$reps==3 & rowCdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1)
  ,xlab="",ylab="",xaxt='n',yaxt='n')
axis(4, las=1)
lines(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$reps==3 & rowCdf$dB==1], col='green1',lwd=5*4)
points(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$reps==3 & rowCdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowCdf$Z[rowCdf$reps==3 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$reps==3 & rowCdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$reps==8 & rowCdf$dB==1], col='green4',lwd=2*4)
lines(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$eFDR_Q[rowCdf$reps==8 & rowCdf$dB==1], col='green4',lwd=2*4)
points(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$reps==8 & rowCdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowCdf$Z[rowCdf$reps==8 & rowCdf$dB==1], y=rowCdf$eFDR_T[rowCdf$reps==8 & rowCdf$dB==1], col='royalblue4',lwd=2*4)
abline(0,1, col='red')
    legText <- c("QuASAR-MPRA","Reads 1/2","Reads 1/5","Reads 1/10","T-test","Reads 1/2","Reads 1/5","Reads 1/10")
    # create a blank graph -- automatically scales -1 to +1 on both axes 
    barplot(0,0, axes=FALSE) 
    legend(x=0.8, y=0, 
      title=paste0("M=60 dB=1"),
           legend=legText, 
           col=c(NA,"green1","green3","green4",NA,"royalblue1","royalblue3","royalblue4"),
           lwd=c(NA,5*4,3*4,2*4,NA,5*4,3*4,2*4),
           xjust=1, 
           yjust=0.5) 
plot(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$power_Q[rowDdf$reads.div==2 & rowDdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1),
  ylab="",yaxt='n')
mtext("p.adj threshold",side=1,line=2.5)
lines(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$power_Q[rowDdf$reads.div==2 & rowDdf$dB==1], col='green1',lwd=5*4)
points(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$power_T[rowDdf$reads.div==2 & rowDdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$power_T[rowDdf$reads.div==2 & rowDdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$power_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$power_Q[rowDdf$reads.div==10 & rowDdf$dB==1], col='green4',lwd=2*4)
lines(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$power_Q[rowDdf$reads.div==10 & rowDdf$dB==1], col='green4',lwd=2*4)
points(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$power_T[rowDdf$reads.div==10 & rowDdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$power_T[rowDdf$reads.div==10 & rowDdf$dB==1], col='royalblue4',lwd=2*4)
plot(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$eFDR_Q[rowDdf$reads.div==2 & rowDdf$dB==1], col='green1',lwd=5*4,xlim=c(0,1), ylim=c(0, 1),
  ylab="",yaxt='n')
axis(4, las=1)
mtext("p.adj threshold",side=1,line=2.5)
lines(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$eFDR_Q[rowDdf$reads.div==2 & rowDdf$dB==1], col='green1',lwd=5*4)
points(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$eFDR_T[rowDdf$reads.div==2 & rowDdf$dB==1], col='royalblue1',lwd=5*4)
lines(x=rowDdf$Z[rowDdf$reads.div==2 & rowDdf$dB==1], y=rowDdf$eFDR_T[rowDdf$reads.div==2 & rowDdf$dB==1], col='royalblue1',lwd=5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_Q[rowAdf$M==60 & rowAdf$dB==1], col='green3',lwd=3.5*4)
points(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
lines(x=rowAdf$Z[rowAdf$M==60 & rowAdf$dB==1], y=rowAdf$eFDR_T[rowAdf$M==60 & rowAdf$dB==1], col='royalblue3',lwd=3.5*4)
points(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$eFDR_Q[rowDdf$reads.div==10 & rowDdf$dB==1], col='green4',lwd=2*4)
lines(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$eFDR_Q[rowDdf$reads.div==10 & rowDdf$dB==1], col='green4',lwd=2*4)
points(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$eFDR_T[rowDdf$reads.div==10 & rowDdf$dB==1], col='royalblue4',lwd=2*4)
lines(x=rowDdf$Z[rowDdf$reads.div==10 & rowDdf$dB==1], y=rowDdf$eFDR_T[rowDdf$reads.div==10 & rowDdf$dB==1], col='royalblue4',lwd=2*4)
abline(0,1, col='red')
dev.off()

