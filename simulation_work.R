library(reshape2)
library(QuASAR)
library(qqman)
#/nfs/rprscratch/wwwShare/cindy/
simSimple1 <- function(prop=0.5,M=80,N=100){
	alpha=prop*M
	beta=M-alpha
	theta=rbeta(length(prop),alpha,beta)
	R=rbinom(length(prop),N,theta)
	cbind(R=R,A=N-R)
}
myttest <- function(simRep_null2_dna){
	rna_ref <- log2(simRep_null2_dna[,c(1,3,5,7,9)]/simRep_null2_dna[,"DNA_R"])
	rna_alt <- log2(simRep_null2_dna[,c(2,4,6,8,10)]/simRep_null2_dna[,"DNA_A"])
	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}

#sim_ttest = function(x=rna_sim) {
#	x_R <- transform(x, 
#	log2_ratio1=log2(X1/DNA_R), 
#	log2_ratio2=log2(X3/DNA_R), 
#	log2_ratio3=log2(X5/DNA_R), 
#	log2_ratio4=log2(X7/DNA_R), 
#  	log2_ratio5=log2(X9/DNA_R))
#	x_R$ID <- 1:nrow(x_R)
#	x_A <- transform(x,
#  	log2_ratio1=log2(X2/DNA_A), 
#  	log2_ratio2=log2(X4/DNA_A), 
#	log2_ratio3=log2(X6/DNA_A), 
#	log2_ratio4=log2(X8/DNA_A), 
# 	log2_ratio5=log2(X10/DNA_A))
#	x_A$ID <- 1:nrow(x_A)
#	x_R_melt <- melt(x_R, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#	names(x_R_melt)[3] <- "log2R_ratio"
#	x_A_melt <- melt(x_A, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#	names(x_A_melt)[3] <- "log2A_ratio"
#	x_melt <- merge(x_R_melt, x_A_melt, by=c("ID","variable"), all=T)
#	x_ttest <- ddply(x_melt, c("ID"), function(y) {
#  	pval_skew = tryCatch(
#    t.test(y$log2R_ratio,y$log2A_ratio, paired=TRUE)$p.value, error=function(y) NA) 
#  	}) 
#	names(x_ttest)[c(2)] <- c("pval_ttest")
#	x_ttest$padj_ttest <- p.adjust(x_ttest$pval_ttest, method="BH")
#	cbind(x,x_ttest)
#}

#DNA
#load("~/piquelab/cindy/bitStarr/Tewhey/summed.RNA.LCL1_ttest.RData")
load("C:/Users/cakal/Google Drive/Lab/Tewhey/revisions/summed.RNA.LCL1_ttest.RData")
dd <- unique(mpra_Filt[,c(1,13:18)])
DNA_sim <- simSimple1(prop=dd$DNA_prop,M=200,N=10000)
DNA_prop <- DNA_sim[,"R"]/(DNA_sim[,"R"]+DNA_sim[,"A"])
DNA_simm <- cbind(matrix(unlist(DNA_sim[,"R"]),ncol=1),matrix(unlist(DNA_sim[,"A"]),ncol=1))
DNA_sim_df <- data.frame(DNA_simm)
colnames(DNA_sim_df) <- c("DNA_R","DNA_A")
simRep_null <- replicate(5,simSimple1(prop=DNA_prop,M=100,N=round((dd$R+dd$A)/5+1)))
dimnames(simRep_null)[[3]] <- paste0("rep",1:5)
simRep_null2 <- matrix(unlist(simRep_null), ncol = 10, byrow = FALSE)
simRep_null2 <- data.frame(simRep_null2)
simRep_null2$R_sum <- rowSums(simRep_null2[,c(1,3,5,7,9)])
simRep_null2$A_sum <- rowSums(simRep_null2[,c(2,4,6,8,10)])
simRep_null2_dna <- cbind(simRep_null2, DNA_sim_df,DNA_prop)
fitQ_null <- fitQuasarMpra(simRep_null2_dna$R_sum,simRep_null2_dna$A_sum,simRep_null2_dna$DNA_prop)

myttest <- function(simRep_null2_dna){
	rna_ref <- log2(simRep_null2_dna[,c(1,3,5,7,9)]/simRep_null2_dna[,"DNA_R"])
	rna_alt <- log2(simRep_null2_dna[,c(2,4,6,8,10)]/simRep_null2_dna[,"DNA_A"])
	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}

myttest_trial <- function(simRep_null2_dna){
	rna_ref <- log2(simRep_null2_dna[ ,grepl("R", names(simRep_null2_dna))]/simRep_null2_dna[,"DNA_R"])
	rna_alt <- log2(simRep_null2_dna[ ,grepl("A", names(simRep_null2_dna))]/simRep_null2_dna[,"DNA_A"])
	apply(rna_ref-rna_alt,1,function(y){tryCatch(t.test(y)$p.value, error=function(x) NA)})
}

pval_ttest <- myttest(simRep_null2_dna)
z <- cbind(simRep_null2_dna,pval_ttest)
z$padj_ttest <- p.adjust(z$pval_ttest, method="BH")

pdf("qqplot_sim_LCL1_null_quasarM100_new.pdf")
qq(fitQ_null$pval3,col ="green")
y <- sort(-log10( z$pval_ttest))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')
legend("topleft", legend=c("QuASAR", "ttest"), 
    fill=c("green", "blue"))
dev.off()

rna_prop1 <- simRep_null2$R_sum/(simRep_null2$A_sum+simRep_null2$R_sum)
pdf("rnaprop_lcl1.pdf")
hist(rna_prop1)
dev.off()

#ind2
load("~/piquelab/cindy/bitStarr/Tewhey/summed.RNA.LCL2_ttest.RData")
dd <- unique(mpra_2_Filt[,-c(2:7,10:11)])
DNA_sim <- simSimple1(prop=dd$DNA_prop,M=200,N=10000)
DNA_prop <- DNA_sim[,"R"]/(DNA_sim[,"R"]+DNA_sim[,"A"])
DNA_simm <- cbind(matrix(unlist(DNA_sim[,"R"]),ncol=1),matrix(unlist(DNA_sim[,"A"]),ncol=1))
DNA_sim_df <- data.frame(DNA_simm)
colnames(DNA_sim_df) <- c("DNA_R","DNA_A")
simRep_null <- replicate(5,simSimple1(prop=DNA_prop,M=100,N=round((dd$R+dd$A)/3+1)))
dimnames(simRep_null)[[3]] <- paste0("rep",1:5)
simRep_null2 <- matrix(unlist(simRep_null), ncol = 10, byrow = FALSE)
simRep_null2 <- data.frame(simRep_null2)
simRep_null2$R_sum <- rowSums(simRep_null2[,c(1,3,5,7,9)])
simRep_null2$A_sum <- rowSums(simRep_null2[,c(2,4,6,8,10)])
simRep_null2_dna <- cbind(simRep_null2, DNA_sim_df,DNA_prop)
fitQ_null <- fitQuasarMpra(simRep_null2_dna$R_sum,simRep_null2_dna$A_sum,simRep_null2_dna$DNA_prop)

simttest_null <- sim_ttest(x=simRep_null2_dna)


pdf("qqplot_sim_LCL2_null_quasarM100_new.pdf")
qq(fitQ_null$pval3,col ="green")
y <- sort(-log10( simttest_null$pval_ttest))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')
legend("topleft", legend=c("QuASAR", "ttest"), 
    fill=c("green", "blue"))
dev.off()

rna_prop2 <- simRep_null2$R_sum/(simRep_null2$A_sum+simRep_null2$R_sum)
pdf("rnaprop_lcl2.pdf")
hist(rna_prop2)
dev.off()

#############reps 
#if want to sim varying counts
#RNA
MPRA_counts <- read.table(file=gzfile("~/piquelab/cindy/bitStarr/Tewhey/GSE75661_79k_collapsed_counts.txt.gz"), sep='\t', header=T)
MPRA_counts$Oligo <- as.character(MPRA_counts$Oligo)
alt <- MPRA_counts[grep('alt', MPRA_counts$Oligo),]
RC <- MPRA_counts[ grep('RC', MPRA_counts$Oligo), ]
R <- MPRA_counts[grep('A',MPRA_counts$Oligo),]
mpra <- transform(MPRA_counts, Allele_class=ifelse(Oligo %in% R$Oligo, "R", "A"), alt_hap=ifelse(Oligo %in% alt$Oligo, "1", "0"), direction =ifelse(Oligo %in% RC$Oligo, "-", "+"))
x <- strsplit(mpra$Oligo, "_")
mpra$rsID <- sapply(x, function(y) { y[1] })
mpra <- transform(mpra, ID=paste(rsID,alt_hap,direction, sep="_"))

rna <- mpra[,-c(2:6,12:19)]
rna_counts_fw <- subset(rna, direction=="+")
lcl1 <- rna_counts_fw[,-c(3:6)]
lcl2 <- rna_counts_fw[,-c(2,4:6)]
lcl3 <- rna_counts_fw[,-c(2:3,5:6)]
lcl4 <- rna_counts_fw[,-c(2:4,6)]
lcl5 <- rna_counts_fw[,-c(3:5)]

dna <- mpra[,-c(7:19)]
dna <- transform(dna, dna=Plasmid_r1+Plasmid_r2+Plasmid_r3+Plasmid_r4+Plasmid_r5)
mpra_fw <- subset(dna, direction=="+")
mpra_counts <- dcast(mpra_fw, ID ~ Allele_class,value.var="dna",fun.aggregate=sum) 
mpra_counts <- transform(mpra_counts, DNA_prop=R/(R+A))
dna_prop <- mpra_counts$DNA_prop
DNA_sim <- simSimple1(prop=dna_prop,M=200,N=10000)
DNA_simm <- cbind(matrix(unlist(DNA_sim[,"R"]),ncol=1),matrix(unlist(DNA_sim[,"A"]),ncol=1))
DNA_sim_df <- data.frame(DNA_simm)
colnames(DNA_sim_df) <- c("DNA_R","DNA_A")
DNA_simF <- subset(DNA_sim_df, DNA_R>=5 & DNA_A>=5)
DNA_simF <- subset(DNA_simF, DNA_R+DNA_A>=100)
DNA_prop <- DNA_simF$DNA_R/(DNA_simF$DNA_R+DNA_simF$DNA_A)
pdf("LCL1_dnaprop_nullsim_filt.pdf")
hist(DNA_prop)
dev.off()

lcl1 <- dcast(lcl1, ID ~ Allele_class, value.var="NA12878_r1",fun.aggregate=sum)
lcl2 <- dcast(lcl2, ID ~ Allele_class, value.var="NA12878_r2",fun.aggregate=sum)
lcl3 <- dcast(lcl3, ID ~ Allele_class, value.var="NA12878_r3",fun.aggregate=sum)
lcl4 <- dcast(lcl4, ID ~ Allele_class, value.var="NA12878_r4",fun.aggregate=sum)
lcl5 <- dcast(lcl5, ID ~ Allele_class, value.var="NA12878_r5",fun.aggregate=sum)

lcl1 <- lcl1[1:length(DNA_prop),]
lcl2 <- lcl2[1:length(DNA_prop),]
lcl3 <- lcl3[1:length(DNA_prop),]
lcl4 <- lcl4[1:length(DNA_prop),]
lcl5 <- lcl5[1:length(DNA_prop),]

simRep_null1 <- simSimple1(prop=DNA_prop,M=100,N=(lcl1$R+lcl1$A)+1)
simRep_null2 <- simSimple1(prop=DNA_prop,M=100,N=(lcl2$R+lcl2$A)+1)
simRep_null3 <- simSimple1(prop=DNA_prop,M=100,N=(lcl3$R+lcl3$A)+1)
simRep_null4 <- simSimple1(prop=DNA_prop,M=100,N=(lcl4$R+lcl4$A)+1)
simRep_null5 <- simSimple1(prop=DNA_prop,M=100,N=(lcl5$R+lcl5$A)+1)
simRep_null_comb <- cbind(simRep_null1, simRep_null2, simRep_null3,simRep_null4,simRep_null5)
df <- matrix(unlist(simRep_null_comb), ncol = 10, byrow = FALSE)
df2 <- data.frame(df)
df2$R_sum <- rowSums(df2[,c(1,3,5,7,9)])
df2$A_sum <- rowSums(df2[,c(2,4,6,8,10)])
df2F <- subset(df2, R_sum>=5 & A_sum>=5)
#df2F <- df2F[1:length(DNA_prop),]

fitQ_null_reps <- fitQuasarMpra(df2F$R_sum,df2F$A_sum,DNA_prop)

simttest_null <- sim_ttest(x=df2F,z=DNA_simF)

simRep_null2_dna <- cbind(simRep_null2, DNA_sim_df,DNA_prop)
fitQ_null <- fitQuasarMpra(simRep_null2_dna$R_sum,simRep_null2_dna$A_sum,simRep_null2_dna$DNA_prop)

simttest_null <- sim_ttest(x=simRep_null2_dna)

pdf("qqplot_sim_LCL1_nullreps_M100_fw.pdf")
qq(fitQ_null_reps$pval3,col ="green")
y <- sort(-log10( simttest_null$pval_ttest))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='blue')
legend("topleft", legend=c("QuASAR", "ttest"), 
    fill=c("green", "blue"))
dev.off()

pdf("/nfs/rprscratch/wwwShare/cindy/qqplot_sim_LCL1_nullreps_quasarM10.pdf")
qq(fitQ_null_reps$pval3)
dev.off()

#######################

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

simRepAse20 <- replicate(5,simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/5+1)))
dimnames(simRepAse20)[[3]] <- paste0("rep",1:5)
#simRepAse20m <- matrix(unlist(simRepAse20), ncol = 10, byrow = FALSE)
simRepAse20df <- data.frame(simRepAse20)
colnames(simRepAse20df) <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
simRepAse20df$R_sum <- rowSums(simRepAse20df[,c(1,3,5,7,9)])
simRepAse20df$A_sum <- rowSums(simRepAse20df[,c(2,4,6,8,10)])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

#fitQ <- fitQuasarMpra(rowSums(simRepAse20[,"R",]),rowSums(simRepAse20[,"A",]),DNA_prop)
fitQ <- fitQuasarMpra(simRepAse20df_dna$R_sum,simRepAse20df_dna$A_sum,simRepAse20df_dna$DNA_prop)
#simttest_ase20 <- sim_ttest(x=simRepAse20df_dna)
#df2F <- subset(df2, R_sum>=5 & A_sum>=5)

table(abs(deltaBeta)>0, fitQ$padj_quasar<0.1)
db <- cbind(simRepAse20df_dna, deltaBeta)
#db <- subset(simRepAse20df_dna_db, R_sum>=5 & A_sum>=5)
table(abs(db$deltaBeta)>0, fitQ$padj_quasar<0.1)

pval_ttest <- myttest(simRepAse20df_dna)
z <- cbind(simRep_null2_dna,pval_ttest)
z$padj_ttest <- p.adjust(z$pval_ttest, method="BH")

all <- cbind(db, fitQ, z)
#all$pval_binom <- sapply(1:nrow(all), function(ii){ 
#	pval_binom=binom.test(c(all[ii,"R_sum"],all[ii,"A_sum"]),all[ii,"DNA_prop"])$p.value
#	})
#fisher test
#n <- transform(all, DNA_R=DNA_R+1, DNA_A=DNA_A+1, R=R_sum+1, A=A_sum+1)#, ID=paste(rsID,alt_hap,sep="_"))
##rownames(n) = make.names(n[,19], unique=TRUE)
#n <- n[,c(11:14)]
#pval_fisher <- apply(n,1, function(x) fisher.test(matrix(x,nr=2))$p.value)
##test <- data.frame(rsID=names(test), pval_fisher=test, row.names=NULL)
##test2 <- transform(test, rsID=gsub("\\.", "\\:", rsID))
##x <- strsplit(test2$rsID, "_")
##test2$rsID <- sapply(x, function(y) { y[1] })
##test2$alt_hap <- sapply(x, function(y) { y[2] })
#all2 <- cbind(all,pval_fisher)
#all2$padj_fisher <- p.adjust(all2$pval_fisher,method="BH")
#
##all_original <- cbind(db, fitQ, simttest_ase20)
#all <- all2
#pdf("qqplot_sim_LCL1_ASE20M100_db2_prop_0.001.pdf")
qq(all$pval3+1E-10,col ="white")
y=-log10((all$pval3+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='green', cex=all$deltaBeta[o]*2+0.01)
y=-log10((all$pval_ttest+1E-10))
o <- order(-y)
x <- (-log10( ppoints( length(o))))
points(x=x, y=y[o], pch=20, col='blue', cex=all$deltaBeta[o]*2+0.01)
#y=-log10((all$pval_binom+1E-10))
#o <- order(-y)
#x <- (-log10( ppoints( length(o))))
#points(x=x, y=y[o], pch=20, col='orange', cex=all$deltaBeta[o]*2+0.01)
#y=-log10((all$pval_fisher+1E-10))
#o <- order(-y)
#x <- (-log10( ppoints( length(o))))
#points(x=x, y=y[o], pch=20, col='lightblue', cex=all$deltaBeta[o]*2+0.01)
y <- sort(-log10( all$pval3[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='red' )
y <- sort(-log10( all$pval_ttest[deltaBeta==0]))
x <- sort(-log10( ppoints( length(y))))
points(x=x, y=y, pch=20, col='purple')
#y <- sort(-log10( all$pval_binom[deltaBeta==0]))
#x <- sort(-log10( ppoints( length(y))))
#points(x=x, y=y, pch=20, col='yellow')
#y <- sort(-log10( all$pval_fisher[deltaBeta==0]))
#x <- sort(-log10( ppoints( length(y))))
#points(x=x, y=y, pch=20, col='darkgreen')
#legend("topleft", legend=c("QuASAR", "ttest","binom","QuASAR null", "ttest null","binom null"),
    #fill=c("green", "blue","orange","red","purple","yellow"))
#dev.off()

#ttest
#simRepAse20_R <- transform(simRepAse20, log2_ratio1=log2(simRepAse20[,"R","rep1"]/DNA_sim[,"R"]), log2_ratio2=log2(simRepAse20[,"R","rep2"]/DNA_sim[,"R"]), 
#	log2_raio3=log2(simRepAse20[,"R","rep3"]/DNA_sim[,"R"]), log2_ratio4=log2(simRepAse20[,"R","rep4"]/DNA_sim[,"R"]), 
#  log2_ratio5=log2(simRepAse20[,"R","rep5"]/DNA_sim[,"R"]))
#simRepAse20_R$ID <- 1:nrow(simRepAse20)
#simRepAse20_A <- transform(simRepAse20,
#  log2_ratio1=log2(simRepAse20[,"A","rep1"]/DNA_sim[,"A"]), log2_ratio2=log2(simRepAse20[,"A","rep2"]/DNA_sim[,"A"]), 
#	log2_ratio3=log2(simRepAse20[,"A","rep3"]/DNA_sim[,"A"]), log2_ratio4=log2(simRepAse20[,"A","rep4"]/DNA_sim[,"A"]), 
#  log2_ratio5=log2(simRepAse20[,"A","rep5"]/DNA_sim[,"A"]))
#simRepAse20_A$ID <- 1:nrow(simRepAse20)
#simRepAse20_R_melt <- melt(simRepAse20_R, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#names(simRepAse20_R_melt)[3] <- "log2R_ratio"
#simRepAse20_A_melt <- melt(simRepAse20_A, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#names(simRepAse20_A_melt)[3] <- "log2A_ratio"
#simRepAse20_melt <- merge(simRepAse20_R_melt, simRepAse20_A_melt, by=c("ID","variable"), all=T)
#
##ttest for each snp, try catch to give NAs when not enough coverage? gave error before I added it
#simRepAse20_ttest <- ddply(simRepAse20_melt, c("ID"), function(x) {
#  pval_skew = tryCatch(
#    t.test(x$log2R_ratio,x$log2A_ratio, paired=TRUE)$p.value, error=function(x) NA) 
#  }) 
#names(simRepAse20_ttest)[c(2)] <- c("pval_ttest")
#simRepAse20_ttest$padj_ttest <- p.adjust(simRepAse20_ttest$pval_ttest, method="BH")
#table(abs(deltaBeta)>0, simRepAse20_ttest$padj_ttest<0.5)
#
#        FALSE  TRUE
#  FALSE 16011  1319
#  TRUE      0  1920
#
#
#pclip = function(p){
#  p + min(p[p > 0],na.rm=T)/10;
#}
#
#
#pdf("/nfs/rprscratch/wwwShare/cindy/qqplot_sim_LCL1_ASE20_quasar.pdf")
#qq(fitQ$pval3)
#dev.off()
#
#pdf("/nfs/rprscratch/wwwShare/cindy/qqplot_sim_LCL1_ASE20_Null_quasar.pdf")
#qq(fitQ$pval3[abs(deltaBeta)==0])
#dev.off()
#
#
#
#
#pdf("qqplot_sim_LCL1_ASE20M100_new_nofxn.pdf")
#qq(fitQ_null_reps$pval3,col ="green")
#y <- sort(-log10( simRepAse20_ttest$pval_ttest))
#x <- sort(-log10( ppoints( length(y))))
#points(x=x, y=y, pch=20, col='blue')
#legend("topleft", legend=c("QuASAR", "ttest"),
#    fill=c("green", "blue"))
#dev.off()
#






#https://github.com/ttriche/ssizeRNA-1/blob/master/R/sim.counts.R

#library(edgeR)
#library(MASS)
#sim = function(m,x) {
#	mu <- dna_prop
#	d <- DGEList(m)
#	d <- calcNormFactors(d)
#	d <- estimateCommonDisp(d)
#	d <- estimateTagwiseDisp(d)
#	M <- 1/d$tagwise.dispersion             ## dispersion for each gene
#	alpha <- mu * M
#	beta <- M - alpha
#	N=x
#	theda <- rbeta(length(mu), alpha, beta, ncp = 0)
#	theda2 <- rbinom(length(mu), N, theda)
#    }
## Fitting the QuASAR model
#rna.res <- fitQuasarMpra(rna_cast$R,rna_cast$A,rna_cast$DNA_prop)
#sum(rna.res$padj_quasar<0.1)
##0
#rna_full <- cbind(rna_cast, df$de, rna.res)
#names(rna_full)[5] <- "de"
#pdf("rna_sim_quasar_041317.pdf")
#qq(rna_full$pval3)
#y <- sort(-log10( rna_full$pval3[ rna_full$de==-1 | rna_full$de==1] ))
#x <- sort(-log10( ppoints( length(y))))
#points(x=x, y=y, pch=20, col='red')
#y <- sort(-log10( rna_full$pval3[ rna_full$de==0] ))
#x <- sort(-log10( ppoints( length(y))))
#points(x=x, y=y, pch=20, col='blue')
#dev.off()
#
#save(rna_full, sc, mu_rna, disp_rna, file="sim_rna_041317.RData")
#simData <- simSimple1(dna_prop,M=80,N=100)



#mu <- mu_dna
#M <- disp_dna
#alpha <- mu * M
#beta <- M - alpha
#N=100
##rbeta(n, shape1, shape2, ncp = 0)
#theda <- rbeta(length(mu), alpha, beta, ncp = 0)
##rbinom(n, size, prob)
#theda2 <- rbinom(length(mu), N, theda)
#simReplicates <- function(prop=0.5,M=80,rna_counts){
#	sapply(1:ncol(rna_counts),function(ii){
#		simSimple1(prop,M,rna_counts[,ii])
#		})
#}
#

#MPRA_counts <- read.table(file=gzfile("~/piquelab/cindy/bitStarr/Tewhey/GSE75661_79k_collapsed_counts.txt.gz"), sep='\t', header=T)
#mpra <- MPRA_counts[,c(1,7:11)]
#rna_counts <- ddply(mpra, c("Oligo","alt_hap"), summarize,
#	NA12878_r1=sum(NA12878_r1))
#counts <- mpra[,2:6]
#
#mpra <- MPRA_counts[,c(1:6)]
#dna_counts <- mpra[,2:6]

#geom.mean = function(row){
#  row[row == 0] = 0.1
#  if(length(row) != 0) return(exp(sum(log(row))/length(row)))
#  else return(0)
#}

#mu_dna <- apply(dna_counts, 1, geom.mean)        ## geometric mean for each gene
#d <- DGEList(dna_counts)
#d <- calcNormFactors(d)
#d <- estimateCommonDisp(d)
#d <- estimateTagwiseDisp(d)
#disp_dna <- d$tagwise.dispersion             ## dispersion for each gene
#
#mu_rna <- apply(counts, 1, geom.mean)        ## geometric mean for each gene
#d <- DGEList(counts)
#d <- calcNormFactors(d)
#d <- estimateCommonDisp(d)
#d <- estimateTagwiseDisp(d)
#disp_rna <- d$tagwise.dispersion             ## dispersion for each gene

#arg=list(nGenes=10000,pi0=0.95,group=rep(c(1,2),each=5))
#source("simulation.R")
#sc <- sim.counts(arg,mu=mu_rna,disp = disp_rna,logfc = 1)

#df <- data.frame(lapply(data.frame(t(sapply(sc, `[`))), unlist))
##names <- rownames(df)
#snp<-rep(c(1:10000),each=10)
#df$ID <- paste("snp",snp,sep="")
#library(plyr)
#rna <- ddply(df, c("ID","group"), summarize,
#	counts=sum(counts)
#	)

#rna_other <- ddply(df, c("ID"), summarize,
#	lambda0=mean(lambda0),
#	phi0=mean(phi0),
#	de=mean(de),
#	delta=mean(delta)
#	)
##ttest
#df2_R <- cbind(df2F, DNA_R=DNA_simF[,1])
#df2_R <- transform(df2_R, 
#	log2_ratio1=log2(X1/DNA_R), 
#	log2_ratio2=log2(X3/DNA_R), 
#	log2_ratio3=log2(X5/DNA_R), 
#	log2_ratio4=log2(X7/DNA_R), 
#	log2_ratio5=log2(X9/DNA_R))
#df2_R$ID <- 1:nrow(df2_R)
#df2_A <- cbind(df2F, DNA_A=DNA_simF[,2])
#df2_A <- transform(df2_A,
#  	log2_ratio1=log2(X2/DNA_A), 
#  	log2_ratio2=log2(X4/DNA_A), 
#	log2_ratio3=log2(X6/DNA_A), 
#	log2_ratio4=log2(X8/DNA_A), 
# 	log2_ratio5=log2(X10/DNA_A))
#df2_A$ID <- 1:nrow(df2_A)
#df2_R_melt <- melt(df2_R, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#names(df2_R_melt)[3] <- "log2R_ratio"
#df2_A_melt <- melt(df2_A, id.vars=c("ID"), measure.vars = c("log2_ratio1","log2_ratio2","log2_ratio3","log2_ratio4","log2_ratio5"))
#names(df2_A_melt)[3] <- "log2A_ratio"
#df2_melt <- merge(df2_R_melt, df2_A_melt, by=c("ID","variable"), all=T)
#
##ttest for each snp, try catch to give NAs when not enough coverage? gave error before I added it
#df2_ttest <- ddply(df2_melt, c("ID"), function(x) {
#  pval_skew = tryCatch(
#    t.test(x$log2R_ratio,x$log2A_ratio, paired=TRUE)$p.value, error=function(x) NA) 
#  }) 
#names(df2_ttest)[c(2)] <- c("pval_ttest")
#df2_ttest$padj_ttest <- p.adjust(df2_ttest$pval_ttest, method="BH")
