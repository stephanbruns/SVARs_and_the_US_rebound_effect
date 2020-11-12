rm(list=ls())

#------------------------------------
# This code produces the results for 
#	dcov and ngml
#	3 variable SVAR
#-----------------------------------

#-----------------------------------
# Packages
#-----------------------------------
require(vars)
library(svars)

library(foreach)
library(doSNOW)
cl<-makeCluster(3) #cores used for parallel computing
registerDoSNOW(cl)

library(normtest)

#-----------------------------------
# Global settings
#-----------------------------------
nbs <- 1000 #bootstrap iterations 
hori_Q <- 80 #horizon quarterly data
hori_M <- 80 #horizon monthly data
set.seed(1234) #for exact reproducibility of the results

loadpath <- "... \\Data\\" #change this
savepath <- "... \\Results\\" #change this

#-----------------------------------
# Functions
#-----------------------------------

# These functions explore the robustness of the identified B matrices
UNIB.dc <- function(VAR.est, iter) {

	SVAR.all <- foreach(i = 1:iter) %dopar% {
		
		library(svars)
		id.dc(VAR.est) #Input: Object of vec2var
	}
	
	B.all <- list()
	for (ii in 1:iter) {
		
		B.all[[ii]] <- SVAR.all[[ii]]$B
	}
	
	return(list(uni = unique(B.all), all = SVAR.all))
	
}	


UNIB.ngml <- function(VAR.est, iter) {

	SVAR.all <- foreach(i = 1:iter) %dopar% {
		
		library(svars)
		id.ngml(VAR.est) #Input: Object of vec2var
	}
	
	B.all <- list()
	for (ii in 1:iter) {
		
		B.all[[ii]] <- SVAR.all[[ii]]$B
	}
	
	return(list(uni = unique(B.all), all = SVAR.all))
	
}	

#######---------------------------------------
####### Analysis of quarterly data
#######---------------------------------------
# Please scroll down for the analysis of monthly data

#-----------------------------------
# Data
#-----------------------------------	

dat <- read.csv(paste(loadpath, "Quarterly Deseasonalized Rebound Data.csv", sep=""), sep=";")

Time.Q <- dat$Time #this is needed to analyze subperiods

dat <- dat[,c("Energy", "RealGDP", "RealPE", "IP", "Quality")]

dat <- log(dat)
colnames(dat) <- c("e","y","p","i","q")
dat <- dat*100 # for easier presentation of the B matrix


#-------
# 5 variables VAR with e, y, p, i and q 
#-------
D <- dat

VAR.est <- VAR(D, p = 3)

# Test for normality
jb.norm.test(VAR.est$varresult$e$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$y$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$p$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$i$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$q$residuals, nrepl=2000)

### dcov
# Uniqueness of B
iterations <- 100
unib.dc <- UNIB.dc(VAR.est, iterations)
unib.dc[[1]]

# B
iden <- unib.dc[[2]][[1]]
iden$B[,1] <- iden$B[,1] * (-1)
iden$B[,5] <- iden$B[,5] * (-1)
#iden$B <- iden$B[, c(1,5,2,3,4)]
round(iden$B, 2) #Mixing matrix of dcov in Table 6 

# Bootstrap
boots <- wild.boot(iden, design = "fixed", distr = "rademacher", nboot = nbs, n.ahead = hori_Q)
res <- list("uniB.dc_Q.v5.p3.levels", iden, boots)

save(res, file=paste(savepath, "irf.dc_Q.v5.p3.levels.iter1000.RDATA", sep=""))


### ngml
# Uniqueness of B
iterations <- 100
unib.ngml <- UNIB.ngml(VAR.est, iterations)
unib.ngml[[1]]

# B
iden <- unib.ngml[[2]][[1]]
iden$B[,1] <- iden$B[,1] * (-1)
round(iden$B, 2) #Mixing matrix of ngml in Table 6 

# Bootstrap
boots <- wild.boot(iden, design = "fixed", distr = "rademacher", nboot = nbs, n.ahead = hori_Q)
res <- list("uniB.ngml_Q.v5.p3.levels", iden, boots)

save(res, file=paste(savepath, "irf.ngml_Q.v5.p3.levels.iter1000.RDATA", sep=""))




#######---------------------------------------
####### Analysis of monthly data
#######---------------------------------------

#-----------------------------------
# Data
#-----------------------------------	
dat <- read.csv(paste(loadpath, "Monthly Deseasonalized Rebound Data.csv", sep=""), sep=";")

dat <- data.frame(dat[,c("Energy", "RealGDP", "RealPE", "IP", "Quality", "HDD", "CDD", "Population")])
colnames(dat) <- c("e","y","p","i","q","h","c", "pop")


dat <- log(dat)
dat <- dat*100

#----------------------------------------
# 5 variables VAR with e, y, p, i and q 
#----------------------------------------
D <- dat[,c("e","y","p","i","q")]

VAR.est <- VAR(D, p = 3)

# Test for normality
jb.norm.test(VAR.est$varresult$e$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$y$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$p$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$i$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$q$residuals, nrepl=2000)

### dcov
# Uniqueness of B
iterations <- 100
unib.dc <- UNIB.dc(VAR.est, iterations)
unib.dc[[1]]

# B
iden <- unib.dc[[2]][[1]]
iden$B[,1] <- iden$B[,1] * (-1)
round(iden$B, 2) #Mixing matrix of dcov in Table 5 

# Bootstrap
boots <-  wild.boot(iden, design = "fixed", distr = "rademacher", nboot = nbs, n.ahead = hori_M)
res <- list("uniB.dc_M.v5.p3.levels", iden, boots)

save(res, file=paste(savepath, "irf.dc_M.v5.p3.levels.iter1000.RDATA", sep=""))


###ngml
# Uniqueness of B
iterations <- 100
unib.ngml <- UNIB.ngml(VAR.est, iterations)
unib.ngml[[1]]

# B
iden <- unib.ngml[[2]][[1]]
iden$B[,1] <- iden$B[,1] * (-1)
round(iden$B, 2) #Mixing matrix of ngml in Table 5 

#Bootstrap
boots <- wild.boot(iden, design = "fixed", distr = "rademacher", nboot = nbs, n.ahead = hori_M)
res <- list("uniB.ngml_M.v5.p3.levels", iden, boots)

save(res, file=paste(savepath, "irf.ngml_M.v5.p3.levels.iter1000.RDATA", sep=""))




#######---------------------------------------
####### Estimation of rebound effects for both monthly and quarterly data
#######---------------------------------------

#-------------------------------------------------------------------------
# Rebound function - calculates the rebound effect for 1, 2, 4 and 6 years
#-------------------------------------------------------------------------

REBOUND <- function(res, nbs, hori) {

	#hori1
	re1 <- ((res[[3]][[1]]$irf[hori[1],2] - res[[3]][[1]]$irf[1,2]) / abs(res[[3]][[1]]$irf[1,2]))
	rebs <- NULL
	for (ii in 1:nbs) {
		rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[hori[1],2] - res[[3]][[2]][[ii]]$irf[1,2]) / abs(res[[3]][[2]][[ii]]$irf[1,2]))
	}
	req1 <- quantile(rebs, prob=c(0.05,0.95))

	#hori2
	re2 <- ((res[[3]][[1]]$irf[hori[2],2] - res[[3]][[1]]$irf[1,2]) / abs(res[[3]][[1]]$irf[1,2]))
	rebs <- NULL
	for (ii in 1:nbs) {
		rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[hori[2],2] - res[[3]][[2]][[ii]]$irf[1,2]) / abs(res[[3]][[2]][[ii]]$irf[1,2]))
	}
	req2 <- quantile(rebs, prob=c(0.05,0.95))
	
	#hori3
	re3 <- ((res[[3]][[1]]$irf[hori[3],2] - res[[3]][[1]]$irf[1,2]) / abs(res[[3]][[1]]$irf[1,2]))
	rebs <- NULL
	for (ii in 1:nbs) {
		rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[hori[3],2] - res[[3]][[2]][[ii]]$irf[1,2]) / abs(res[[3]][[2]][[ii]]$irf[1,2]))
	}
	req3 <- quantile(rebs, prob=c(0.05,0.95))
	
	#hori4
	re4 <- ((res[[3]][[1]]$irf[hori[4],2] - res[[3]][[1]]$irf[1,2]) / abs(res[[3]][[1]]$irf[1,2]))
	rebs <- NULL
	for (ii in 1:nbs) {
		rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[hori[4],2] - res[[3]][[2]][[ii]]$irf[1,2]) / abs(res[[3]][[2]][[ii]]$irf[1,2]))
	}
	req4 <- quantile(rebs, prob=c(0.05,0.95))

	return(list(hori1 = c(re1, req1), hori2 = c(re2, req2), hori3 = c(re3, req3), hori4 = c(re4, req4)))	
}

#-------
# dcov
#-------

# Create vector with filenames for loop
filename <- c("irf.dc_Q.v5.p3.levels.iter1000.RDATA",
				"irf.dc_M.v5.p3.levels.iter1000.RDATA")
				
# Create table to collect rebound effects
reb_table <- data.frame(matrix(ncol=9, nrow=length(filename)))
colnames(reb_table) <- c("Name", "1 year","0.90 CI", "2 years","0.90 CI", "4 years", "0.90 CI", "6 years", "0.90 CI")

# Create list to collect B matrices
B_list <- list()

# Loop for the different models
for (i in 1:length(filename)) {
	load(paste(savepath, filename[i], sep=""))

	# 1, 2, 4 and 6 years
	if (i < 2) {horiz <- c(4,8,16,24)}else{horiz <- c(12,24,48,72)} 
	
	REB <- REBOUND(res, nbs=1000, hori=horiz) 
	reb_table[i,1] <- filename[i]
	reb_table[i,2] <- round(REB[[1]][1],2)
	reb_table[i,3] <- paste("[", round(REB[[1]][2],2),",",round(REB[[1]][3],2),"]", sep="")	
	
	reb_table[i,4] <- round(REB[[2]][1],2)
	reb_table[i,5] <- paste("[", round(REB[[2]][2],2),",",round(REB[[2]][3],2),"]", sep="")	

	reb_table[i,6] <- round(REB[[3]][1],2)
	reb_table[i,7] <- paste("[", round(REB[[3]][2],2),",",round(REB[[3]][3],2),"]", sep="")	
	
	reb_table[i,8] <- round(REB[[4]][1],2)
	reb_table[i,9] <- paste("[", round(REB[[4]][2],2),",",round(REB[[4]][3],2),"]", sep="")	
		
	# IRF
	pdf(paste(savepath, filename[i],".pdf", sep=""))
	print(plot(res[[3]]))
	dev.off()

	# B matrix
	B_list[[i]] <- list(filename[i], res[[2]]$B)

}

# Rebound effects for dcov (Table 4)
reb_table

	
#-------
# ngml
#-------

# Create vector with filenames for loop
filename <- c("irf.ngml_Q.v5.p3.levels.iter1000.RDATA",
				"irf.ngml_M.v5.p3.levels.iter1000.RDATA")
				
# Create table to collect rebound effects
reb_table <- data.frame(matrix(ncol=9, nrow=length(filename)))
colnames(reb_table) <- c("Name", "1 year","0.90 CI", "2 years","0.90 CI", "4 years", "0.90 CI", "6 years", "0.90 CI")

# Create table to collect B matrices
B_list <- list()

# Loop for the different models
for (i in 1:length(filename)) {
	load(paste(savepath, filename[i], sep=""))

	# 1, 2, 4 and 6 years
	if (i < 2) {horiz <- c(4,8,16,24)}else{horiz <- c(12,24,48,72)} 
	
	REB <- REBOUND(res, nbs=1000, hori=horiz) 
	reb_table[i,1] <- filename[i]
	reb_table[i,2] <- round(REB[[1]][1],2)
	reb_table[i,3] <- paste("[", round(REB[[1]][2],2),",",round(REB[[1]][3],2),"]", sep="")	
	
	reb_table[i,4] <- round(REB[[2]][1],2)
	reb_table[i,5] <- paste("[", round(REB[[2]][2],2),",",round(REB[[2]][3],2),"]", sep="")	

	reb_table[i,6] <- round(REB[[3]][1],2)
	reb_table[i,7] <- paste("[", round(REB[[3]][2],2),",",round(REB[[3]][3],2),"]", sep="")	
	
	reb_table[i,8] <- round(REB[[4]][1],2)
	reb_table[i,9] <- paste("[", round(REB[[4]][2],2),",",round(REB[[4]][3],2),"]", sep="")	
		
	# IRF
	pdf(paste(savepath, filename[i],".pdf", sep=""))
	print(plot(res[[3]]))
	dev.off()

	# B matrix
	B_list[[i]] <- list(filename[i], res[[2]]$B)

}

# Rebound effects for ngml (Table 4)
reb_table



























