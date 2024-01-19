library("ggplot2")
load("pH0_0.8L\\pH0_0.8VL1.Rdata")
sPowerList <- truePosList
sAlfaList  <- falsePosList
sValueList <- valueList
sEffortList<- effortList
sRetryList <- rrList
sFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

finalVal1  <- matrix(finalVal,ncol=1)
finalEff1  <- matrix(finalEff,ncol=1)
finalRR1   <- matrix(finalRR,ncol=1)
finalFP1   <- colSums(finalcFP)
finalTP1   <- colSums(finalcTP)
finalFPP1  <- finalFP1/(finalFP1+finalTP1)

fppDF               <- data.frame(fpp=finalFPP1)
fppDF$SimCondition  <- "1. One Hypothesis"

histDF              <- data.frame(effort=finalEff1,rr=finalRR1)
histDF$SimCondition <- "1. One Hypothesis"

load("pH0_0.8L\\pH0_0.8VL2.Rdata")

rPowerList <- truePosList
rAlfaList  <- falsePosList
rValueList <- valueList
rEffortList<- effortList
rRetryList <- rrList
rFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

finalVal05  <- matrix(finalVal,ncol=1)
finalEff05  <- matrix(finalEff,ncol=1)
finalRR05   <- matrix(finalRR,ncol=1)
finalFP05   <- colSums(finalcFP)
finalTP05   <- colSums(finalcTP)
finalFPP05  <- finalFP05/(finalFP05+finalTP05)

hist2DF     <- data.frame(effort=finalEff05,rr=finalRR05)
hist2DF$SimCondition <- "2. Many Hypotheses, 0.05"
fpp2DF      <- data.frame(fpp=finalFPP05)
fpp2DF$SimCondition  <- "2. Many Hypotheses, 0.05"

load("pH0_0.8L\\pH0_0.8VL3.Rdata")
pPowerList <- truePosList
pAlfaList  <- falsePosList
pValueList <- valueList
pEffortList<- effortList
pRetryList <- rrList
pFPList    <- falsePosList[,2]/(truePosList[,2]+falsePosList[,2])

yrs        <- seq(1,nSimSteps)

finalVal005  <- matrix(finalVal,ncol=1)
finalEff005  <- matrix(finalEff,ncol=1)
finalRR005   <- matrix(finalRR,ncol=1)
finalFP005   <-colSums(finalcFP)
finalTP005   <-colSums(finalcTP)
finalFPP005  <- finalFP005/(finalFP005+finalTP005)

hist3DF     <- data.frame(effort=finalEff005,rr=finalRR005)
hist3DF$SimCondition <- "3. Many Hypotheses, 0.005"
fpp3DF      <- data.frame(fpp=finalFPP005)
fpp3DF$SimCondition  <- "3. Many Hypotheses, 0.005"

histDF <- rbind(histDF,hist2DF,hist3DF)
fppDF <- rbind(fppDF,fpp2DF,fpp3DF)

dev.new(width=12,height=8)
par(mfrow=c(2,2))
plot(yrs,sEffortList[,2],xlab="steps",ylab="Effort",ylim=c(0,100))
points(yrs,rEffortList[,2],col="blue")
points(yrs,pEffortList[,2],col="red")

plot(yrs,sValueList[,2],xlab="steps",ylab="Value",ylim=c(0,400))
points(yrs,rValueList[,2],col="blue")
points(yrs,pValueList[,2],col="red")

plot(yrs,sRetryList[,2],ylim=c(0,40),xlab="steps",ylab="Number of hypotheses")
points(yrs,rRetryList[,2],col="blue")
points(yrs,pRetryList[,2],col="red")
legend(x="topleft",legend=c("1 hypothesis","many hypotheses","p<0.005, many hypotheses"),lty=c(1,1,1),col=c("black","blue","red"),lwd=3)

plot(yrs,sFPList,ylim=c(0,1),xlab="steps",ylab="Fraction of *False Pos* Pubs")
points(yrs,rFPList,col="blue")
points(yrs,pFPList,col="red")

dev.new(width=12,height=8)
p1<-ggplot(histDF, aes(x=effort, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.75, binwidth = 10, boundary = 0)  +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02")) +
xlim(c(0,150))+xlab("Effort")
print(p1)
dev.new(width=12,height=8)
p2<-ggplot(fppDF, aes(x=fpp, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.75, binwidth = 0.025, boundary = 0) +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02")) +
xlim(c(0,1)) + xlab("Fraction of False Positive Publications")
print(p2)
dev.new(width=12,height=8)
p3<-ggplot(histDF, aes(x=rr, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.75, binwidth = 2, boundary = 0) +
scale_color_manual(values=c("#000000", "#090de0","#de0a02" ))+scale_fill_manual(values=c("#404040", "#090de0", "#de0a02"))+
xlim(c(0,40)) + xlab("Number of Hypotheses")
print(p3)

histDFD <- histDF[which(histDF$SimCondition!="1. One Hypothesis"),]
dev.new(width=12,height=8)
p4<-ggplot(histDFD, aes(x=rr, color=SimCondition, fill=SimCondition)) + geom_histogram(position = "identity", alpha = 0.75, binwidth = 2, boundary = 0) +
scale_color_manual(values=c("#090de0","#de0a02" ))+scale_fill_manual(values=c("#090de0", "#de0a02"))+
xlim(c(0,40)) + xlab("Number of Hypotheses")
print(p4)