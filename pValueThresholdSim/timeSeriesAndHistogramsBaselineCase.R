#
# pick a data set
load("pH0_0.8L\\pH0_0.8VL2.Rdata")
#
# orient data for plotting
#
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
#
# generate plots
#
dev.new()
plot(truePosList[,2],ylim=c(0,1),xlab="time step",ylab="True Discovery Rate")
lines(truePosList[,1],col="red")
lines(truePosList[,3],col="red")
dev.new()
plot(falsePosList[,2],ylim=c(0,1),xlab="time step",ylab="False Discovery Rate")
lines(falsePosList[,1],col="red")
lines(falsePosList[,3],col="red")
dev.new()
plot(effortList[,2],xlab="time step",ylab="Effort",ylim=c(0,maxEff))
lines(effortList[,1],col="red")
lines(effortList[,3],col="red")
dev.new()
plot(valueList[,2],xlab="time step",ylab="Value",ylim=c(0,max(valueList)))
lines(valueList[,1],col="red")
lines(valueList[,3],col="red")
dev.new()
plot(rrList[,2],ylim=c(0,maxHypotheses),xlab="time step",ylab="number of hypotheses")
lines(rrList[,1],col="red")
lines(rrList[,3],col="red")

dev.new()
hist(finalEff1,breaks=seq(5,500,by=5),main="Effort histogram",xlab="Effort",ylab="Frequency",xlim=c(0,100))

dev.new()
hist(finalVal1,breaks=seq(0,4500,by=20),main="Value histogram",xlab="Value",ylab="Frequency",xlim=c(0,800))

dev.new()
hist(finalRR1,breaks=seq(0,80,by=1),main="Number of hypotheses histogram",xlab="Number of hypotheses",ylab="Frequency",xlim=c(0,40))

dev.new()
hist(finalFPP1,breaks=seq(0,1,by=0.025),main="False positive rate histogram",xlab="False positive rate",ylab="Frequency",xlim=c(0,1))
