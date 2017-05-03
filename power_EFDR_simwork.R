load("C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/5_rep_varyManddB.RData")

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
	NumTrue <- sum(abs(i$deltaBeta)>0)
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
zz <- lapply(trial,t)

#loop over list containing power and eFDR calculations and plot QuASAR vs T-test
zz[] <- lapply(zz,function(x) replace(x,is.na(x),0))
#plot z by power for QuASAR and TTest over M_dB combinations
lapply(zz, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,".png"))
	plot(x=data$Z, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$power_Q, col='green')
	points(x=data$Z, y=data$power_T, col='blue')
	lines(x=data$Z, y=data$power_T, col='blue')
	dev.off()
	})

#plot eFDR by power for QuASAR and TTest over M_dB combinations
lapply(zz, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDRxPower.png"))
	plot(x=data$eFDR_Q, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="eFDR",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$eFDR_Q, y=data$power_Q, col='green')
	points(x=data$eFDR_T, y=data$power_T, col='blue')
	lines(x=data$eFDR_T, y=data$power_T, col='blue')
	dev.off()
	})

#needs work!
#plot eFDR by power separating z into line style
dat <- do.call("rbind", zz)
zzm <- matrix(unlist(dat),ncol=7,byrow = FALSE)
zzdf <- data.frame(zzm)
colnames(zzdf) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")
zzdf_sub <- subset(zzdf, M==10)
M_label=zzdf_sub$M[1]
#dB_label=zzdf_sub$dB[1]
png(paste0("M",M_label,"eFDRxPower_combineddB.png"))
plot(x=zzdf$eFDR_Q[zzdf$M==10], y=zzdf$power_Q[zzdf$M==10], col='green',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="eFDR",ylab="power",main=paste0("M",M_label," ","zzdf$dB=",dB_label))
lines(x=zzdf$eFDR_Q[zzdf$M==10], y=zzdf$power_Q[zzdf$M==10], col='green',lty=3)
points(x=zzdf$eFDR_T[zzdf$M==10], y=zzdf$power_T[zzdf$M==10], col='blue',lty=3)
lines(x=zzdf$eFDR_T[zzdf$M==10], y=zzdf$power_T[zzdf$M==10], col='blue',lty=3)
points(x=zzdf$eFDR_Q[zzdf$M==60], y=zzdf$power_Q[zzdf$M==60], col='green',lty=1)
lines(x=zzdf$eFDR_Q[zzdf$M==60], y=zzdf$power_Q[zzdf$M==60], col='green',lty=1)
points(x=zzdf$eFDR_T[zzdf$M==60], y=zzdf$power_T[zzdf$M==60], col='blue',lty=1)
lines(x=zzdf$eFDR_T[zzdf$M==60], y=zzdf$power_T[zzdf$M==60], col='blue',lty=1)
points(x=zzdf$eFDR_Q[zzdf$M==100], y=zzdf$power_Q[zzdf$M==100], col='green',lty=2)
lines(x=zzdf$eFDR_Q[zzdf$M==100], y=zzdf$power_Q[zzdf$M==100], col='green',lty=2)
points(x=zzdf$eFDR_T[zzdf$M==100], y=zzdf$power_T[zzdf$M==100], col='blue',lty=2)
lines(x=zzdf$eFDR_T[zzdf$M==100], y=zzdf$power_T[zzdf$M==100], col='blue',lty=2)
#points(x=zzdf$eFDR_Q[zzdf$dB==2], y=zzdf$power_Q[zzdf$dB==2], col='green',,lty=3)
#lines(x=zzdf$eFDR_Q[zzdf$dB==2], y=zzdf$power_Q[zzdf$dB==2], col='green',lty=3)
#points(x=zzdf$eFDR_T[zzdf$dB==2], y=zzdf$power_T[zzdf$dB==2], col='blue',lty=3)
#lines(x=zzdf$eFDR_T[zzdf$dB==2], y=zzdf$power_T[zzdf$dB==2], col='blue',lty=3)
dev.off()