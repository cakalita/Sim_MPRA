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

M10_db2
M10_db1
M10_db0_5
M10_db0_1
M60_db2
M60_db1
M60_db0_5
M60_db0_1
M100_db2
M100_db1
M100_db0_5
M100_db0_1

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

#data <- data.frame()
#for (i in c("M10_db2","M10_db1","M10_db0_5","M10_db0_1","M60_db2","M60_db1","M60_db0_5","M60_db0_1","M100_db2","M100_db1","M100_db0_5","M100_db0_1")) {
#q <- table(abs(deltaBeta[i])>0, padj_quasar[i]<0.1)
#t <- table(abs(deltaBeta[i])>0, padj_ttest[i]<0.1)
##M <- substr(i,2,nchar(i)-2)
#z <- data.frame(power_Q=q[2,2]/(q[2,2]+q[2,1]), eFDR_Q=q[1,2]/(q[2,2]+q[1,2]),
#	power_T=t[2,2]/(t[2,2]+t[2,1]), eFDR_T=t[1,2]/(t[2,2]+t[1,2]))
#data <- rbind(data,cbind(z,group=i))
#}
list1 <- list(M10_db2,M10_db1,M10_db0_5,M10_db0_1,M60_db2,M60_db1,M60_db0_5,M60_db0_1,M100_db2,M100_db1,M100_db0_5,M100_db0_1)

#trial <-sapply(list1, function(i) {
#	q <- table(abs(i$deltaBeta)>0, i$padj_quasar<0.1,exclude=NULL)
#	t <- table(abs(i$deltaBeta)>0, i$padj_ttest<0.1,exclude=NULL)
#	z <- data.frame(
#		power_Q=q[2,2]/(q[2,2]+q[2,1]), 
#		eFDR_Q=q[1,2]/(q[2,2]+q[1,2]),
#		power_T=t[2,2]/(t[2,2]+t[2,1]), 
#		eFDR_T=t[1,2]/(t[2,2]+t[1,2]))
#	}, simplify=TRUE)

#repeat this for however many y (padj thresholds)
trial_old <-sapply(list1, function(i,y=0.1) {
	q <- table(abs(i$deltaBeta)>0, i$padj_quasar<y,exclude=NULL)
	t <- table(abs(i$deltaBeta)>0, i$padj_ttest<y,exclude=NULL)
	z <- data.frame(
		power_Q=q[2,2]/(q[2,2]+q[2,1]), 
		eFDR_Q=q[1,2]/(q[2,2]+q[1,2]),
		power_T=t[2,2]/(t[2,2]+t[2,1]), 
		eFDR_T=t[1,2]/(t[2,2]+t[1,2]))
	}, simplify=TRUE)


zvec <- c(0.01,0.02,0.05,0.08,(1:10)/10);

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
trial1 <- lapply(trial,t)
# library(stringr)
#z2 <- lapply(z, function(x) str_trim(x))
#newm <- matrix(unlist(trial_old), ncol = 4, byrow = FALSE)
#newdf <- data.frame(newm)
#row.names(newdf) <- c("M10_db2","M10_db1","M10_db0_5","M10_db0_1","M60_db2","M60_db1","M60_db0_5","M60_db0_1","M100_db2","M100_db1","M100_db0_5","M100_db0_1")
#colnames(newdf) <- c("power_Q", "eFDR_Q", "power_T", "eFDR_T")
#newdf_0_1 <- transform(newdf, M=c(rep(10,4),rep(60,4),rep(100,4)),dB=rep(c(2,1,0.5,0.1),3),Z=0.1)
#
#z <- list(newdf_0,newdf_0_1,newdf_0_2,newdf_0_3,newdf_0_4,newdf_0_5,newdf_0_6,newdf_0_7,newdf_0_8,newdf_0_9,newdf_1)
zz <- list(lapply(z[1], subset, M==10 & dB==0.1),
	lapply(z[2], subset,  M==10 & dB==0.1),
	lapply(z[3], subset,  M==10 & dB==0.1),
	lapply(z[4], subset,  M==10 & dB==0.1),
	lapply(z[5], subset,  M==10 & dB==0.1),
	lapply(z[6], subset,  M==10 & dB==0.1),
	lapply(z[7], subset,  M==10 & dB==0.1),
	lapply(z[8], subset,  M==10 & dB==0.1),
	lapply(z[9], subset,  M==10 & dB==0.1),
	lapply(z[10], subset, M==10 & dB==0.1),
	lapply(z[11], subset, M==10 & dB==0.1))
zz <- lapply(trial1, subset, M==10 & dB==2)

zz_unlist <- data.frame(lapply(data.frame(t(sapply(zz, `[`))), unlist))
zz_unlist_t <- t(zz_unlist)
# col 7 is the padj threshold (Z), column 1 is power_Q, col 3 is power_T
plot(x=zz_unlist_t[,7], y=zz_unlist_t[,1], col='green',ylim=c(0, 1),xlab="padj",ylab="power",main="M10_db0.1")
lines(x=zz_unlist_t[,7], y=zz_unlist_t[,1], col='green')
points(x=zz_unlist_t[,7], y=zz_unlist_t[,3], col='blue')
lines(x=zz_unlist_t[,7], y=zz_unlist_t[,3], col='blue')

zz[] <- lapply(trial1,function(x) replace(x,is.na(x),0))
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


  q <- table(abs(M60_db0_1$deltaBeta)>0, M60_db0_1$padj_quasar<0.1,exclude=NULL)
  t <- table(abs(M60_db0_1$deltaBeta)>0, M60_db0_1$padj_ttest<0.1,exclude=NULL)

q[2,2]/(q[2,2]+q[2,1]) 

t[2,2]/(t[2,2]+t[2,1]) 
