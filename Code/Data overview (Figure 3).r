rm(list=ls())

#---------------------------------------------
# This code produces the data graph (Figure 3)
#---------------------------------------------

#-----------------------------------
# Paths
#-----------------------------------	
loadpath <- "... \\Data\\" #change this
savepath <- "... \\Results\\" #change this

#---------------------------------------
# Data
#---------------------------------------

# quarterly data
dat.Q <- read.csv(paste(loadpath, "Quarterly Deseasonalized Rebound Data.csv", sep=""), sep=";")
Time.Q <- dat.Q$Time
Time.Q <- as.character(Time.Q)
Time.Q.axis <- c(as.character(Time.Q), "2016-4", "2017-1") #axis
dat.Q <- dat.Q[,c("Energy", "RealGDP", "RealPE", "IP", "Quality")]
dat.Q <- log(dat.Q)

# monthly data
dat.M <- read.csv(paste(loadpath, "Monthly Deseasonalized Rebound Data.csv", sep=""), sep=";")
Time.M <- dat.M$Time
Time.M <- as.character(Time.M)
Time.M.axis <- c(as.character(Time.M), "2016-11", "2016-12", "2017-1") # axis
dat.M <- dat.M[,c("Energy", "RealGDP", "RealPE", "IP", "Quality")]
dat.M <- log(dat.M)

#---------------------------------------
# Figure 3
#---------------------------------------
pdf(paste(savepath, "Data_figure.pdf", sep=""), width=10, height=12)
par(mfrow=c(nrows=5, ncols=2))

#Energy
matplot(1:length(Time.M), dat.M$Energy, t="l", ylab="Energy", xlab="", axes=F, main="Monthly")
axis(2)
axis(side=1,at=seq(1, length(Time.M.axis), 12),labels=Time.M.axis[seq(1, length(Time.M.axis), 12)], xlab=F)

matplot(1:length(Time.Q), dat.Q$Energy, t="l", ylab="Energy", xlab="", axes=F, main="Quarterly")
axis(2)
axis(side=1,at=seq(1, length(Time.Q.axis), 4),labels=Time.Q.axis[seq(1, length(Time.Q.axis), 4)], xlab=F)

#GDP
matplot(1:length(Time.M), dat.M$RealGDP, t="l", ylab="GDP", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.M.axis), 12),labels=Time.M.axis[seq(1, length(Time.M.axis), 12)], xlab=F)

matplot(1:length(Time.Q), dat.Q$RealGDP, t="l", ylab="GDP", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.Q.axis), 4),labels=Time.Q.axis[seq(1, length(Time.Q.axis), 4)], xlab=F)

#Energy price
matplot(1:length(Time.M), dat.M$RealPE, t="l", ylab="Energy price", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.M.axis), 12),labels=Time.M.axis[seq(1, length(Time.M.axis), 12)], xlab=F)

matplot(1:length(Time.Q), dat.Q$RealPE, t="l", ylab="Energy price", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.Q.axis), 4),labels=Time.Q.axis[seq(1, length(Time.Q.axis), 4)], xlab=F)

#IP
matplot(1:length(Time.M), dat.M$IP, t="l", ylab="Industrial production", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.M.axis), 12),labels=Time.M.axis[seq(1, length(Time.M.axis), 12)], xlab=F)

matplot(1:length(Time.Q), dat.Q$IP, t="l", ylab="Industrial production", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.Q.axis), 4),labels=Time.Q.axis[seq(1, length(Time.Q.axis), 4)], xlab=F)

#Quality
matplot(1:length(Time.M), dat.M$Quality, t="l", ylab="Energy quality", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.M.axis), 12),labels=Time.M.axis[seq(1, length(Time.M.axis), 12)], xlab=F)

matplot(1:length(Time.Q), dat.Q$Quality, t="l", ylab="Energy quality", xlab="", axes=F)
axis(2)
axis(side=1,at=seq(1, length(Time.Q.axis), 4),labels=Time.Q.axis[seq(1, length(Time.Q.axis), 4)], xlab=F)

dev.off()




