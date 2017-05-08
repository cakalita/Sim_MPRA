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

#plot z by eFDR for QuASAR and TTest over M_dB combinations
lapply(zz, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDR.png"))
	plot(x=data$Z, y=data$eFDR_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="eFDR",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$eFDR_Q, col='green')
	points(x=data$Z, y=data$eFDR_T, col='blue')
	lines(x=data$Z, y=data$eFDR_T, col='blue')
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
png(paste0("M",M_label,"dB1_zxPower_combineddB.png"))
plot(x=zzdf$Z[zzdf$M==10 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==10 & zzdf$dB==0.5], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="power",main=paste0("M=varied dB=0.5"))
lines(x=zzdf$Z[zzdf$M==10 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==10 & zzdf$dB==0.5], col='lightgreen',lty=3)
points(x=zzdf$Z[zzdf$M==10 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==10 & zzdf$dB==0.5], col='lightblue',lty=3)
lines(x=zzdf$Z[zzdf$M==10 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==10 & zzdf$dB==0.5], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==0.5], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==0.5], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==0.5], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==0.5], col='blue',lty=1)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==100 & zzdf$dB==0.5], col='darkgreen',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==100 & zzdf$dB==0.5], col='darkgreen',lty=2)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==100 & zzdf$dB==0.5], col='darkblue',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==100 & zzdf$dB==0.5], col='darkblue',lty=2)
#points(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$power_Q[zzdf$dB==2 & zzdf$dB==1], col='green',,lty=3)
#lines(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$power_Q[zzdf$dB==2 & zzdf$dB==1], col='green',lty=3)
#points(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$power_T[zzdf$dB==2 & zzdf$dB==1], col='blue',lty=3)
#lines(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$power_T[zzdf$dB==2 & zzdf$dB==1], col='blue',lty=3)
dev.off()
plot(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==0.5], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="power",main=paste0("M=60 dB=varied"))
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==0.5], col='lightgreen',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==0.5], col='lightblue',lty=3)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==0.5], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==0.5], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==2], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==2], col='darkgreen',lty=2)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==2], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==2], col='darkgreen',lty=2)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==2], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==2], col='darkblue',lty=2)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==2], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==2], col='darkblue',lty=2)
#eFDR plots
plot(x=zzdf$Z[zzdf$M==10 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==10 & zzdf$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="eFDR",main=paste0("M=varied dB=1"))
lines(x=zzdf$Z[zzdf$M==10 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==10 & zzdf$dB==1], col='lightgreen',lty=3)
points(x=zzdf$Z[zzdf$M==10 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==10 & zzdf$dB==1], col='lightblue',lty=3)
lines(x=zzdf$Z[zzdf$M==10 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==10 & zzdf$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==1], col='darkgreen',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==1], col='darkgreen',lty=2)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==1], col='darkblue',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==1], col='darkblue',lty=2)
abline(0,1, col='red')
#points(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$dB==2 & zzdf$dB==1], col='green',,lty=3)
#lines(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$dB==2 & zzdf$dB==1], col='green',lty=3)
#points(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$dB==2 & zzdf$dB==1], col='blue',lty=3)
#lines(x=zzdf$Z[zzdf$dB==2 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$dB==2 & zzdf$dB==1], col='blue',lty=3)
dev.off()
plot(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==0.5], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="eFDR",main=paste0("M=100 dB=varied"))
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==0.5], col='lightgreen',lty=3)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==0.5], col='lightblue',lty=3)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==0.5], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==0.5], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==1], col='blue',lty=1)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)
lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)
abline(0,1, col='red')

#repeat for ind 2
load("C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/3_rep_varyManddB.RData")

#should be able to remove this paragraph when I go back and save the objects with padj in varyManddB.R
R3_M10_db2 <-  transform(R3_M10_db2, M=10, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M10_db1 <-  transform(R3_M10_db1, M=10, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M10_db0_5 <-  transform(R3_M10_db0_5, M=10, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M10_db0_1 <-  transform(R3_M10_db0_1, M=10, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M60_db2 <-  transform(R3_M60_db2, M=60, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M60_db1 <-  transform(R3_M60_db1, M=60, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M60_db0_5 <-  transform(R3_M60_db0_5, M=60, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M60_db0_1 <-  transform(R3_M60_db0_1, M=60, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M100_db2 <-  transform(R3_M100_db2, M=100, dB=2, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M100_db1 <-  transform(R3_M100_db1, M=100, dB=1, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M100_db0_5 <-  transform(R3_M100_db0_5, M=100, dB=0.5, padj_ttest=p.adjust(pval_ttest, method="BH"))
R3_M100_db0_1 <-  transform(R3_M100_db0_1, M=100, dB=0.1, padj_ttest=p.adjust(pval_ttest, method="BH"))

#make list object of different simulations (changing variance and effect size (similar to fold change))
list2 <- list(R3_M10_db2,R3_M10_db1,R3_M10_db0_5,R3_M10_db0_1,R3_M60_db2,R3_M60_db1,R3_M60_db0_5,R3_M60_db0_1,R3_M100_db2,R3_M100_db1,R3_M100_db0_5,R3_M100_db0_1)

#set up series of padj thresholds
zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);
#loop over set of simulations (possibly could work on this and add in the actual simulation with loop over M and dB conditions more automatically)
trial2 <-lapply(list2, function(i) {
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
zz2 <- lapply(trial2,t)

#loop over list containing power and eFDR calculations and plot QuASAR vs T-test
zz2[] <- lapply(zz2,function(x) replace(x,is.na(x),0))
#plot z by power for QuASAR and TTest over M_dB combinations
lapply(zz2, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"padjxpower_IND2.png"))
	plot(x=data$Z, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$power_Q, col='green')
	points(x=data$Z, y=data$power_T, col='blue')
	lines(x=data$Z, y=data$power_T, col='blue')
	dev.off()
	})

#plot z by eFDR for QuASAR and TTest over M_dB combinations
lapply(zz2, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDR_IND2.png"))
	plot(x=data$Z, y=data$eFDR_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="eFDR",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$eFDR_Q, col='green')
	points(x=data$Z, y=data$eFDR_T, col='blue')
	lines(x=data$Z, y=data$eFDR_T, col='blue')
	dev.off()
	})

#plot eFDR by power for QuASAR and TTest over M_dB combinations
lapply(zz2, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDRxPower_IND2.png"))
	plot(x=data$eFDR_Q, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="eFDR",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$eFDR_Q, y=data$power_Q, col='green')
	points(x=data$eFDR_T, y=data$power_T, col='blue')
	lines(x=data$eFDR_T, y=data$power_T, col='blue')
	dev.off()
	})

#plot ind 1 vs ind 2 dB1 M60
dat2 <- do.call("rbind", zz2)
zzm2 <- matrix(unlist(dat2),ncol=7,byrow = FALSE)
zzdf2 <- data.frame(zzm2)
colnames(zzdf2) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")

plot(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="power",main=paste0("M=60 dB=1 Ind 1 (5 rep) vs Ind 2 (3 rep)"))
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3)
points(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$power_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$power_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
#points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$power_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
#lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$power_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
#points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$power_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)
#lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$power_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)

plot(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="Z",ylab="eFDR",main=paste0("M=60 dB=1 Ind 1 (5 rep) vs Ind 2 (3 rep)"))
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_Q[zzdf2$M==60 & zzdf2$dB==1], col='lightgreen',lty=3)
points(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
lines(x=zzdf2$Z[zzdf2$M==60 & zzdf2$dB==1], y=zzdf2$eFDR_T[zzdf2$M==60 & zzdf2$dB==1], col='lightblue',lty=3)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_Q[zzdf$M==60 & zzdf$dB==1], col='green',lty=1)
points(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
lines(x=zzdf$Z[zzdf$M==60 & zzdf$dB==1], y=zzdf$eFDR_T[zzdf$M==60 & zzdf$dB==1], col='blue',lty=1)
#points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
#lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_Q[zzdf$M==100 & zzdf$dB==2], col='darkgreen',lty=2)
#points(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)
#lines(x=zzdf$Z[zzdf$M==100 & zzdf$dB==2], y=zzdf$eFDR_T[zzdf$M==100 & zzdf$dB==2], col='darkblue',lty=2)
abline(0,1, col='red')

#needs work!
#plot eFDR by power separating z into line style
dat <- do.call("rbind", zz2)
zz2m <- matrix(unlist(dat),ncol=7,byrow = FALSE)
zz2df <- data.frame(zz2m)
colnames(zz2df) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")
zz2df_sub <- subset(zz2df, M==10)
M_label=zz2df_sub$M[1]
#dB_label=zz2df_sub$dB[1]
png(paste0("M",M_label,"eFDRxPower_combineddB_IND2.png"))
plot(x=zz2df$eFDR_Q[zz2df$M==10], y=zz2df$power_Q[zz2df$M==10], col='green',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="eFDR",ylab="power",main=paste0("M",M_label," ","zz2df$dB=",dB_label))
lines(x=zz2df$eFDR_Q[zz2df$M==10], y=zz2df$power_Q[zz2df$M==10], col='green',lty=3)
points(x=zz2df$eFDR_T[zz2df$M==10], y=zz2df$power_T[zz2df$M==10], col='blue',lty=3)
lines(x=zz2df$eFDR_T[zz2df$M==10], y=zz2df$power_T[zz2df$M==10], col='blue',lty=3)
points(x=zz2df$eFDR_Q[zz2df$M==60], y=zz2df$power_Q[zz2df$M==60], col='green',lty=1)
lines(x=zz2df$eFDR_Q[zz2df$M==60], y=zz2df$power_Q[zz2df$M==60], col='green',lty=1)
points(x=zz2df$eFDR_T[zz2df$M==60], y=zz2df$power_T[zz2df$M==60], col='blue',lty=1)
lines(x=zz2df$eFDR_T[zz2df$M==60], y=zz2df$power_T[zz2df$M==60], col='blue',lty=1)
points(x=zz2df$eFDR_Q[zz2df$M==100], y=zz2df$power_Q[zz2df$M==100], col='green',lty=2)
lines(x=zz2df$eFDR_Q[zz2df$M==100], y=zz2df$power_Q[zz2df$M==100], col='green',lty=2)
points(x=zz2df$eFDR_T[zz2df$M==100], y=zz2df$power_T[zz2df$M==100], col='blue',lty=2)
lines(x=zz2df$eFDR_T[zz2df$M==100], y=zz2df$power_T[zz2df$M==100], col='blue',lty=2)
#points(x=zz2df$eFDR_Q[zz2df$dB==2], y=zz2df$power_Q[zz2df$dB==2], col='green',,lty=3)
#lines(x=zz2df$eFDR_Q[zz2df$dB==2], y=zz2df$power_Q[zz2df$dB==2], col='green',lty=3)
#points(x=zz2df$eFDR_T[zz2df$dB==2], y=zz2df$power_T[zz2df$dB==2], col='blue',lty=3)
#lines(x=zz2df$eFDR_T[zz2df$dB==2], y=zz2df$power_T[zz2df$dB==2], col='blue',lty=3)
dev.off()

#repeat for 8 reps
load("C:/Users/Cindy/Google Drive/Lab/Tewhey/revisions/8_rep_varyManddB.RData")

#should be able to remove this paragraph when I go back and save the objects with padj in varyManddB.R 
## think I want to keep it
M10_db2 <-  transform(M10_db2, M=10, dB=2)
M10_db1 <-  transform(M10_db1, M=10, dB=1)
M10_db0_5 <-  transform(M10_db0_5, M=10, dB=0.5)
M10_db0_1 <-  transform(M10_db0_1, M=10, dB=0.1)
M60_db2 <-  transform(M60_db2, M=60, dB=2)
M60_db1 <-  transform(M60_db1, M=60, dB=1)
M60_db0_5 <-  transform(M60_db0_5, M=60, dB=0.5)
M60_db0_1 <-  transform(M60_db0_1, M=60, dB=0.1)
M100_db2 <-  transform(M100_db2, M=100, dB=2)
M100_db1 <-  transform(M100_db1, M=100, dB=1)
M100_db0_5 <-  transform(M100_db0_5, M=100, dB=0.5)
M100_db0_1 <-  transform(M100_db0_1, M=100, dB=0.1)

#make list object of different simulations (changing variance and effect size (similar to fold change))
list8 <- list(M10_db2,M10_db1,M10_db0_5,M10_db0_1,M60_db2,M60_db1,M60_db0_5,M60_db0_1,M100_db2,M100_db1,M100_db0_5,M100_db0_1)

#set up series of padj thresholds
zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);
#loop over set of simulations (possibly could work on this and add in the actual simulation with loop over M and dB conditions more automatically)
trial8 <-lapply(list8, function(i) {
	#i <- list8[[1]]
	NumTrue <- sum(abs(i$deltaBeta)>0, na.rm=T)
	aux <- sapply(zvec,function(z){
		NumDiscQ <- sum(i$padj_quasar<z,na.rm=T)
		NumDiscT <- sum(i$padj_ttest<z,na.rm=T) 
		NumTrueDiscQ <- sum(i$padj_quasar<z & abs(i$deltaBeta)>0, na.rm=T) 
		NumTrueDiscT <- sum(i$padj_ttest<z & abs(i$deltaBeta)>0, na.rm=T) 
		M <- i$M[1]
		deltaBeta = i$dB[1]
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
zz8 <- lapply(trial8,t)

#loop over list containing power and eFDR calculations and plot QuASAR vs T-test
zz8[] <- lapply(zz8,function(x) replace(x,is.na(x),0))
#plot z by power for QuASAR and TTest over M_dB combinations
lapply(zz8, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"padjxpower_8rep.png"))
	plot(x=data$Z, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$power_Q, col='green')
	points(x=data$Z, y=data$power_T, col='blue')
	lines(x=data$Z, y=data$power_T, col='blue')
	dev.off()
	})

#plot z by eFDR for QuASAR and TTest over M_dB combinations
lapply(zz8, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDR_8rep.png"))
	plot(x=data$Z, y=data$eFDR_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="padj",ylab="eFDR",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$Z, y=data$eFDR_Q, col='green')
	points(x=data$Z, y=data$eFDR_T, col='blue')
	lines(x=data$Z, y=data$eFDR_T, col='blue')
	dev.off()
	})

#plot eFDR by power for QuASAR and TTest over M_dB combinations
lapply(zz8, function(x) {
	data=data.frame(x)
	M_label=data$M[1]
	dB_label=data$dB[1]
	png(paste0("M",M_label,"_dB",dB_label,"eFDRxPower_rep8.png"))
	plot(x=data$eFDR_Q, y=data$power_Q, col='green',xlim=c(0,1), ylim=c(0, 1),xlab="eFDR",ylab="power",main=paste0("M",M_label," ","dB=",dB_label))
	lines(x=data$eFDR_Q, y=data$power_Q, col='green')
	points(x=data$eFDR_T, y=data$power_T, col='blue')
	lines(x=data$eFDR_T, y=data$power_T, col='blue')
	dev.off()
	})

#needs work!
#plot eFDR by power separating z into line style
dat <- do.call("rbind", zz8)
zz8m <- matrix(unlist(dat),ncol=7,byrow = FALSE)
zz8df <- data.frame(zz8m)
colnames(zz8df) <- c("power_Q","eFDR_Q" ,"power_T","eFDR_T","Z" ,"M","dB")
zz8df_sub <- subset(zz8df, M==10)
M_label=zz8df_sub$M[1]
#dB_label=zz8df_sub$dB[1]
png(paste0("M",M_label,"eFDRxPower_combineddB_8rep.png"))
plot(x=zz8df$eFDR_Q[zz8df$M==10], y=zz8df$power_Q[zz8df$M==10], col='green',lty=3,xlim=c(0,1), ylim=c(0, 1),
	xlab="eFDR",ylab="power",main=paste0("M",M_label," ","zz8df$dB=",dB_label))
lines(x=zz8df$eFDR_Q[zz8df$M==10], y=zz8df$power_Q[zz8df$M==10], col='green',lty=3)
points(x=zz8df$eFDR_T[zz8df$M==10], y=zz8df$power_T[zz8df$M==10], col='blue',lty=3)
lines(x=zz8df$eFDR_T[zz8df$M==10], y=zz8df$power_T[zz8df$M==10], col='blue',lty=3)
points(x=zz8df$eFDR_Q[zz8df$M==60], y=zz8df$power_Q[zz8df$M==60], col='green',lty=1)
lines(x=zz8df$eFDR_Q[zz8df$M==60], y=zz8df$power_Q[zz8df$M==60], col='green',lty=1)
points(x=zz8df$eFDR_T[zz8df$M==60], y=zz8df$power_T[zz8df$M==60], col='blue',lty=1)
lines(x=zz8df$eFDR_T[zz8df$M==60], y=zz8df$power_T[zz8df$M==60], col='blue',lty=1)
points(x=zz8df$eFDR_Q[zz8df$M==100], y=zz8df$power_Q[zz8df$M==100], col='green',lty=2)
lines(x=zz8df$eFDR_Q[zz8df$M==100], y=zz8df$power_Q[zz8df$M==100], col='green',lty=2)
points(x=zz8df$eFDR_T[zz8df$M==100], y=zz8df$power_T[zz8df$M==100], col='blue',lty=2)
lines(x=zz8df$eFDR_T[zz8df$M==100], y=zz8df$power_T[zz8df$M==100], col='blue',lty=2)
#points(x=zz8df$eFDR_Q[zz8df$dB==2], y=zz8df$power_Q[zz8df$dB==2], col='green',,lty=3)
#lines(x=zz8df$eFDR_Q[zz8df$dB==2], y=zz8df$power_Q[zz8df$dB==2], col='green',lty=3)
#points(x=zz8df$eFDR_T[zz8df$dB==2], y=zz8df$power_T[zz8df$dB==2], col='blue',lty=3)
#lines(x=zz8df$eFDR_T[zz8df$dB==2], y=zz8df$power_T[zz8df$dB==2], col='blue',lty=3)
dev.off()