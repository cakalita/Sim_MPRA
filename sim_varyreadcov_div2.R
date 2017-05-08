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
simRepAse20 <- replicate(5,simSimple1(prop=RNA_prop,M=10,N=round((dd$R+dd$A)/2+1)))
simRepAse20df <- data.frame(simRepAse20)
simRepAse20df$r_sum <- rowSums(simRepAse20df[ ,grepl("R", names(simRepAse20df))])
simRepAse20df$a_sum <- rowSums(simRepAse20df[ ,grepl("A", names(simRepAse20df))])
simRepAse20df_dna <- cbind(simRepAse20df, DNA_sim_df,DNA_prop)

fitQ <- fitQuasarMpra(simRepAse20df_dna$r_sum,simRepAse20df_dna$a_sum,simRepAse20df_dna$DNA_prop)
db <- cbind(simRepAse20df_dna, deltaBeta)
pval_ttest_trial <- data.frame(pval_ttest=myttest_trial(simRepAse20df_dna))
pval_ttest_trial$padj_ttest <- p.adjust(pval_ttest_trial$pval_ttest, method="BH")

M10_db2 <- cbind(db, fitQ, pval_ttest)
#M10_db1 <- cbind(db, fitQ, pval_ttest)
#M10_db0_5 <- cbind(db, fitQ, pval_ttest)
#M60_db2 <- cbind(db, fitQ, pval_ttest)
#M60_db1 <- cbind(db, fitQ, pval_ttest)
#M60_db0_5 <- cbind(db, fitQ, pval_ttest)
#M100_db2 <- cbind(db, fitQ, pval_ttest)
#M100_db1 <- cbind(db, fitQ, pval_ttest)
#M100_db0_5 <- cbind(db, fitQ, pval_ttest)

save(M10_db2,M10_db1,M10_db0_5,M60_db2,M60_db1,M60_db0_5,
	M100_db2,M100_db1,M100_db0_5,
	file="C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/5_rep_div2_varyManddB.RData")

#should be able to remove this paragraph when I go back and save the objects with padj in varyManddB.R
M10_db2 <-  transform(M10_db2, M=10, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
M10_db1 <-  transform(M10_db1, M=10, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
M10_db0_5 <-  transform(M10_db0_5, M=10, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
M10_db0_1 <-  transform(M10_db0_1, M=10, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))
M60_db2 <-  transform(M60_db2, M=60, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
M60_db1 <-  transform(M60_db1, M=60, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
M60_db0_5 <-  transform(M60_db0_5, M=60, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
M60_db0_1 <-  transform(M60_db0_1, M=60, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))
M100_db2 <-  transform(M100_db2, M=100, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
M100_db1 <-  transform(M100_db1, M=100, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
M100_db0_5 <-  transform(M100_db0_5, M=100, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
M100_db0_1 <-  transform(M100_db0_1, M=100, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))

#make list object of different simulations (changing variance and effect size (similar to fold change))
list1 <- list(M10_db2,M10_db1,M10_db0_5,M10_db0_1,M60_db2,M60_db1,M60_db0_5,M60_db0_1,M100_db2,M100_db1,M100_db0_5,M100_db0_1)

#set up series of padj thresholds
zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);
#loop over set of simulations (possibly could work on this and add in the actual simulation with loop over M and dB conditions more automatically)
trial <-lapply(list1, function(i) {
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
zz_div10 <- lapply(trial,t)

#loop over list containing power and eFDR calculations and plot QuASAR vs T-test
zz_div10[] <- lapply(zz_div10,function(x) replace(x,is.na(x),0))


zz_div2


#plot ind 1 vs ind 2 dB1 M60
dat2 <- do.call("rbind", zz_div2)
zzm2 <- matrix(unlist(dat2),ncol=7,byrow = FALSE)
zzdf_div2 <- data.frame(zzm2)
colnames(zzdf_div2) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")

dat2 <- do.call("rbind", zz_div10)
zzm2 <- matrix(unlist(dat2),ncol=7,byrow = FALSE)
zzdf_div10 <- data.frame(zzm2)
colnames(zzdf_div10) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")

plot(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$power_Q[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="power",main=paste0("M=60 dB=1 varied read coverage (reads / 2 | 5 | 10)"))
lines(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$power_Q[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightgreen',lty=3)
points(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$power_T[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightblue',lty=3)
lines(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$power_T[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$power_Q[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkgreen',lty=2)
lines(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$power_Q[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkgreen',lty=2)
points(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$power_T[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkblue',lty=2)
lines(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$power_T[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkblue',lty=2)

plot(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$eFDR_Q[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="eFDR",main=paste0("M=60 dB=1 varied read coverage (reads / 2 | 5 | 10)"))
lines(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$eFDR_Q[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightgreen',lty=3)
points(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$eFDR_T[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightblue',lty=3)
lines(x=zzdf_div2$Z[zzdf_div2$M==60 & zzdf_div2$dB==1], y=zzdf_div2$eFDR_T[zzdf_div2$M==60 & zzdf_div2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$eFDR_Q[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkgreen',lty=2)
lines(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$eFDR_Q[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkgreen',lty=2)
points(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$eFDR_T[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkblue',lty=2)
lines(x=zzdf_div10$Z[zzdf_div10$M==60 & zzdf_div10$dB==1], y=zzdf_div10$eFDR_T[zzdf_div10$M==60 & zzdf_div10$dB==1], col='darkblue',lty=2)
abline(0,1, col='red')

