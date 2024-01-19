rm(list=ls())
library(dplyr)
library("tictoc")
tic()
#
# Specify simulation parameters
#
myFileName         <- "pH0_0.5VL3.Rdata"       # file for saving output
nSimSteps          <- 50000                  # number of simulation time steps
nMC                <- 100                    # number of Monte Carlo replication runs
nLabs              <- 2000                   # number of researchers in the sim

nHypotheses        <- 1                      # initial number of hypotheses (p-hacking model)
alfa               <- 0.005                  # p-value threshold for testing hypotheses  
maxHypotheses      <- 50                     # upper limit on number of hypotheses

pH0True   		   <- 0.5                    # prior probability of true null
deltMin            <- 0.2                    # minimum true effect size for false null

effort             <- 20                     # initial effort
maxEff      	   <- 500                    # maximum allowable effort
minEff             <- 5                      # minimum allowable effort

LM                 <- 0.1                    # fraction of additional labs to draw from in evolution
deathRate          <- 0.002                  # probability of lab death
delayTime          <- 100                    # time step at which evolution begins

sigma_RR           <- 1                      # standard deviation for nHypotheses evolution
sigma_e            <- 10                     # standard deviation for effort evolution

eta                <- 0.01;                  # shape parameter for effort-to-experiment-attempt function

#
# Initialize arrays for the simulation
#
lab_id             <- 1:nLabs
testMat            <- matrix(0,ncol=nLabs,nrow=maxHypotheses)  # testMat is an indicator for hypothesis tests for each lab
deltMat            <- testMat;                                 # deltMat is a matrix for the true effect sizes, labs X nHypotheses
#
# labs is a data frame containing the records for all living researchers
#
labs               <- data.frame(matrix(ncol=11, nrow= nLabs))
colnames(labs)     <- c('age','id','effort','truePos','nHypotheses','value','falsePos','runExp','cumTruePos','cumFalsePos','pub_rate')
#
# save quartiles, across MC runs, for various quantities of interest
#
falsePosList   <- matrix(0,ncol=3,nrow=nSimSteps)
truePosList    <- matrix(0,ncol=3,nrow=nSimSteps)
falseNegList   <- matrix(0,ncol=3,nrow=nSimSteps)
trueNegList    <- matrix(0,ncol=3,nrow=nSimSteps)
effortList     <- matrix(0,ncol=3,nrow=nSimSteps)
valueList      <- matrix(0,ncol=3,nrow=nSimSteps)
rrList         <- matrix(0,ncol=3,nrow=nSimSteps)
mutList        <- matrix(0,ncol=3,nrow=nSimSteps)

#
# save final results for each living lab and each MC run
#
finalRR   <- matrix(0,nrow=nLabs,ncol=nMC)
finalEff  <- matrix(0,nrow=nLabs,ncol=nMC)
finalVal  <- matrix(0,nrow=nLabs,ncol=nMC)
finalFP   <- matrix(0,nrow=nLabs,ncol=nMC)
finalTP   <- matrix(0,nrow=nLabs,ncol=nMC)
finalAge  <- matrix(0,nrow=nLabs,ncol=nMC)
finalDelt <- matrix(0,nrow=nLabs,ncol=nMC)
finalcFP  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcTP  <- matrix(0,nrow=nLabs,ncol=nMC)
finalcExp <- matrix(0,nrow=nLabs,ncol=nMC)
#
#  indxMat is used for indexing the researchers' number of hypotheses over time
#
indxMat    <- rep(seq(1,maxHypotheses),times=nLabs);
indxMat    <- t(matrix(indxMat,ncol=nLabs))
#
# monte carlo replication loop
#
for (iMC in 1:nMC){

#
#  initialize the lab dataframe
#
    print(paste("iMC = ",iMC,sep=""))
	labs$age               <- 0
	labs$id                <- seq(1,nLabs)
	labs$effort            <- effort
	labs$truePos           <- 0
	labs$nHypotheses       <- nHypotheses
	labs$value             <- 0
	labs$falsePos          <- 0
	labs$runExp            <- 0
	labs$cumTruePos        <- 0
	labs$cumFalsePos       <- 0
	
#
#  simulate the base effect sizes for the labs 
#
	delta1                 <- rgamma(nLabs,1,rate = -log(1-pH0True)/deltMin)
	myDelta                <- matrix(rep(delta1,times=maxHypotheses),ncol=maxHypotheses)
#
# time stepping loop
#
	for (i in 1:nSimSteps) {

	  h       <- exp(-(labs$effort-minEff)*eta)  # probability of conducting experiment this time step 
	  w       <- runif(nLabs,0,1)

	  
	  runExp       <- ifelse(h>=w,1,0);          # simulate whether or not to run experiment
	  labs$runExp  <- labs$runExp+runExp;
#
#     Set the t-critical values for alfa and Sidak correction
#	  
	  alfSidak   <- 1-(1-alfa)^(1/labs$nHypotheses)     
	  t_alfSidak <- qt(1-alfSidak/2, 2*labs$effort-2)   #critical T for multiple testing
	  t_alfa     <- qt(1-alfa/2, 2*labs$effort-2)       #critical T for one-test 

#
#  initialize arrays for hypothesis testing
#
	  detct    <- rep(0,nLabs);
	  truth    <- rep(0,nLabs);
	  detCorr  <- rep(0,nLabs);	  
	  rrMat    <- matrix(rep(labs$nHypotheses,times=maxHypotheses),ncol=maxHypotheses);
	  testMat  <- ifelse(indxMat<=rrMat,1,0)	  
#
#    Simulate effect sizes using the lab effect
#	  
	  myDelta4 <- matrix(rgamma(nLabs,1,scale = myDelta),ncol=1)	  
	  myDelta3 <- matrix(rep(myDelta4,times=maxHypotheses),ncol=maxHypotheses);
	  myDelta2 <- (myDelta3)*testMat
#
#     find the true positives: myTruth is an indicator matrix of false nulls
#
	  myTruth  <- ifelse(myDelta2<=deltMin,0,1)
#
#     generate t-statistics for each lab, each hypothesis
#     non-central t's used for true effect sizes
#	  
	  effMat   <- matrix(rep(labs$effort,times=maxHypotheses),ncol=maxHypotheses); # effect sizes
	  lamda    <- sqrt(effMat)*myDelta2/sqrt(2)                                    # noncentrality
	  tVals    <- abs(rt(rrMat,2*effMat-2,lamda)*testMat)                          # t statistics
	  tMax     <- apply(tVals,1,max)	                                           # maximize over hypotheses  
	  tMaxMat  <- matrix(rep(tMax,times=maxHypotheses),ncol=maxHypotheses);        
	  chkMax   <- ifelse(tMaxMat == tVals,1,0);                                    # which hypothesis meets the max
	  detct    <- tMax>t_alfa                                                      # max beats standard T
	  detCorr  <- tMax>t_alfSidak                                                  # max beats Sidak T
	  truth    <- ifelse(rowSums(chkMax*myTruth)>0,1,0)                            # max was at a false null?
	  
	  truePos <- ifelse(((truth==1)&(detCorr==1)&(h>=w)),1,0)                      # false null, beat Sidak, experiment ran
	  falsPos <- ifelse(((truth==0)&(detct==1)&(h>=w)),1,0)                        # true null, beat T, experiment ran
#
#  record experiment outcomes in labs data frame
#	  
	  labs$falsePos      <- falsPos
	  labs$truePos       <- truePos
	  labs$value         <- labs$value + (truePos + falsPos)
	  labs$age           <- labs$age + 1
	  labs$cumTruePos    <- labs$cumTruePos + truePos
	  labs$cumFalsePos   <- labs$cumFalsePos + falsPos 
	  labs$pub_rate      <- labs$value/labs$age
#
# if time is right, perform evolution
#	  
	  if(i>delayTime){
#
#  who dies 
#	  
	  	deaths      <- runif(nLabs,0,1)
		indxDeath   <- which(deaths<=deathRate);
#
#       sort the labs
#	  
	    L           <- length(indxDeath);                                      # number of deaths
	    L2          <- min(floor(L*(1+LM)),nLabs);                             # number of "best" labs to draw from
	    bestLabs    <- sort.int(labs$value,index.return=TRUE,decreasing=TRUE)  # sort labs in terms of value
	    indxRep     <- sample(bestLabs$ix[1:L2],L,replace=FALSE)               # indices of L labs sampled from the best L*(1+LM) 
	                           
	    labs[indxDeath,]       <- labs[indxRep,]   # initialize new labs with data from "best"
#
#  reset parameters not associated with evolution inheritance
#
	    labs$age[indxDeath]            <-1
	    labs$runExp[indxDeath]         <-0
	    labs$value[indxDeath]          <-0
	    labs$cumFalsePos[indxDeath]    <-0
	    labs$cumTruePos[indxDeath]     <-0
	    labs$pub_rate[indxDeath]       <-0
	    labs$id[indxDeath]             <- seq(max(labs$id)+1,max(labs$id)+L)
#
#  perform evolution in terms of effort and nHypotheses
#	  
	    labs$effort[indxDeath] <- round(labs$effort[indxDeath] + rnorm(L,0,sigma_e))
	    labs$effort            <- ifelse(labs$effort>maxEff,maxEff,labs$effort)
	    labs$effort            <- ifelse(labs$effort<minEff,minEff,labs$effort)
	  
	    labs$nHypotheses[indxDeath] <- round(labs$nHypotheses[indxDeath] + rnorm(L,0,sigma_RR))
	    labs$nHypotheses            <- ifelse(labs$nHypotheses>maxHypotheses,maxHypotheses,labs$nHypotheses)
	    labs$nHypotheses            <- ifelse(labs$nHypotheses<1,1,labs$nHypotheses)  
#
#   set effect size of new labs from old labs
#
#	    myDelta[indxDeath,] <- myDelta[indxRep,]
		deltaD                 <- rgamma(L,1,rate = -log(1-pH0True)/deltMin)
		myDelta[indxDeath,]    <- matrix(rep(deltaD,times=maxHypotheses),ncol=maxHypotheses)	    
	  }
	  
	  truePosList[i,]  <- truePosList[i,] + c(sum(labs$truePos)/nLabs-1/sqrt(nLabs),sum(labs$truePos)/nLabs,sum(labs$truePos)/nLabs+1/sqrt(nLabs))/nMC 
	  
	  falsePosList[i,]   <- falsePosList[i,] + c(sum(labs$falsePos)/nLabs-1/sqrt(nLabs),sum(labs$falsePos)/nLabs,sum(labs$falsePos)/nLabs+1/sqrt(nLabs))/nMC  
	  
	  trueNegList[i,]  <- trueNegList[i,] + c(sum(runExp*(1-labs$falsePos))/nLabs-1/sqrt(nLabs),sum(runExp*(1-labs$falsePos))/nLabs,sum(runExp*(1-labs$falsePos))/nLabs+1/sqrt(nLabs))/nMC 
	  
	  falseNegList[i,]   <- falsePosList[i,] + c(sum(runExp*(1-labs$truePos))/nLabs-1/sqrt(nLabs),sum(runExp*(1-labs$truePos))/nLabs,sum(runExp*(1-labs$truePos))/nLabs+1/sqrt(nLabs))/nMC   
	  
	  effortList[i,] <- effortList[i,] + quantile(labs$effort,probs = c(0.25,0.5,0.75))/nMC
	  valueList[i,]  <- valueList[i,] + quantile(labs$value,probs = c(0.25,0.5,0.75))/nMC
	  rrList[i,]     <- rrList[i,] + quantile(labs$nHypotheses,probs = c(0.25,0.5,0.75))/nMC

	}
	medDelt         <- 0*labs$effort
	for (iii in 1:nLabs){
		 medDelt[iii]<- median(myDelta[iii,(1:labs$nHypotheses[iii])])
	}
	finalRR[,iMC]   <- labs$nHypotheses
	finalEff[,iMC]  <- labs$effort
	finalVal[,iMC]  <- labs$value
	finalFP[,iMC]   <- labs$falsePos
	finalTP[,iMC]   <- labs$truePos
	finalAge[,iMC]  <- labs$age
	finalDelt[,iMC] <- medDelt
	finalcFP[,iMC]  <- labs$cumFalsePos
	finalcTP[,iMC]  <- labs$cumTruePos
	finalcExp[,iMC] <- labs$runExp
}

print(labs)

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
plot(rrList[,2],ylim=c(0,maxHypotheses),xlab="time step",ylab="re-try rate")
lines(rrList[,1],col="red")
lines(rrList[,3],col="red")
toc()
save.image(file=myFileName)
