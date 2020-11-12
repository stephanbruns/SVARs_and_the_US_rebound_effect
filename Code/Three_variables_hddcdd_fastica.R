rm(list=ls())


#------------------------------------
# This code produces the results for 
#	fastICA
#	3 variables SVAR (energy-gdp-prices) 
# including two exogenous variables: HDD (heating degree days) and CDD (cooling degree days)
#
#-------------------------------------------------------------------------

# package path, ****** please CHANGE next line accordingly!!! **** :
ppath<-".../Replication package/" #change this

#-----------------------------------
# Packages and source files
#-----------------------------------
require(vars)
require(normtest)
require(fastICA)
library(svars)
library(pbapply)
library(clue)
library(steadyICA)
library(expm)

source(paste(ppath, "Code/sources/sourcedir.R", sep=""))
sourceDir(paste(ppath, "Code/sources/SVARS_mod/", sep=""), FALSE)
sourceDir(paste(ppath, "Code/sources/maincodes/", sep=""), FALSE)

#-----------------------------------
# Global settings
#-----------------------------------
nbs <- 1000 #bootstrap iterations 
hori_Q <- 80 #horizon quarterly data
hori_M <- 80 #horizon monthly data
sseed<-42 #for exact reproducibility of the results

loadpath <- paste(ppath, "Data/", sep="") 
savepath <- paste(ppath, "Results/", sep="") 

# function to estimate fastICA:
fAp<-function(ures){
  set.seed(sseed)
  X<-t(ures)
  icares <- fastICA(t(X), nrow(X),tol=1e-14, maxit=3000, verbose=FALSE) #Doris: tol=1e-14
  W <- t((icares$K) %*% (icares$W)) 
  A <-solve(W) # A is the mixing matrix
  eres<- t(W %*% X) # 
  # we first rescale the mixing matrix so that the independent components (i.e. structural shocks) have variance equal to one
  Ascaled<-A
  eresscaled<-eres
  for (i in 1:nrow(A)){
    eresscaled[,i]<-eres[,i]/sd(eres[,i])
    Ascaled[,i]<-A[,i]*sd(eres[,i])    
  }
  aba<-abs(Ascaled)
  cc1<-which(aba == max(aba), arr.ind = TRUE) # coordinate of the matrix where is the max entry
  aba[cc1[1],]<-0
  aba[,cc1[2]]<-0
  cc2<-which(aba == max(aba), arr.ind = TRUE)
  aba[cc2[1],]<-0
  aba[,cc2[2]]<-0
  cc3<-which(aba == max(aba), arr.ind = TRUE)
  CC<-rbind(cc1,cc2,cc3)
  e_s<-CC[CC[,1]==1,2] #energyshock
  y_s<-CC[CC[,1]==2,2] #gdpshock
  p_s<-CC[CC[,1]==3,2] #priceshock
  #perme<-c(e_s, (1:k)[-e_s])
  perme<-c(e_s, y_s, p_s)
  Ap<-Ascaled[,perme]
  colnames(Ap) <-c("e","y","p")
  if(Ap[1,1]>0){Ap[,1]<- -Ap[,1]}
  for (i in 2:k){if(Ap[i,i]<0){Ap[,i]<- -Ap[,i]}}
  Ap
}

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

#----------------------------------
# 3 variables VAR with e, y and p + 2 exogenous variables: hdd + cdd
#----------------------------------
D <- dat[,c("e","y","p")]
k<-ncol(D)
hdd <- ts(dat[,"h"], start = 1992, frequency = 12)
cdd <- ts(dat[,"c"], start = 1992, frequency = 12)
exovar=cbind(hdd, cdd)

VAR.est <- vars::VAR(D, p = 3, type="const", exogen=exovar) # VARX with exogenous variabes
ures<-resid(VAR.est) # reduced form VAR residuals

BX<-matrix(c(VAR.est$varresult$e$coefficients["hdd"],
             VAR.est$varresult$y$coefficients["hdd"],
             VAR.est$varresult$p$coefficients["hdd"],
             VAR.est$varresult$e$coefficients["cdd"],
             VAR.est$varresult$y$coefficients["cdd"],
             VAR.est$varresult$p$coefficients["cdd"]),
           3,2)

# Test for normality
jb.norm.test(VAR.est$varresult$e$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$y$residuals, nrepl=2000)
jb.norm.test(VAR.est$varresult$p$residuals, nrepl=2000)

### fastICA estimation
iden<-id.fastica(VAR.est, sseed=sseed)
Ap<-fAp(ures)
iden$B<-Ap # this is the mixing matrix
round(iden$B,2) #Mixing matrix of fastICA in Table 1 (confidence intervals are calculated below in the rebound function)


# confidence intervals
#boots <- wild.boot_mod(iden, rademacher = T, nboot = nbs, n.ahead = hori_Q)
boots <- wild.boot_mod_exgvar(iden, rademacher = T, nboot = nbs, n.ahead = hori_M,
                              exog_varb = exovar, Bex=BX)
q10 <- function(x) {quantile(x, probs=c(0.05,0.95))}
contemp_sig <- matrix(ncol = ncol(boots$bootstrap[[1]]$irf), nrow = length(boots$bootstrap))
colnames(contemp_sig) <- colnames(boots$bootstrap[[1]]$irf)
for (j in 1:length(boots$bootstrap)) {
  contemp_sig[j,] <- unlist(boots$bootstrap[[j]]$irf[1,])
}
c.i.<-apply(contemp_sig, 2, q10)
round(c.i.,2) #### Confidence intervals of fastICA in Table 1

#save data
res <- list("uniB.fastica_M.v3.hddcdd.p3.levels", iden, boots)
save(res, file=paste(savepath, "irf.fastica_M.v3.hddcdd.p3.levels.iter1000.RDATA", sep=""))

#######---------------------------------------
####### Estimation of rebound effects 
#######---------------------------------------


#-------------------------------------------------------------------------
# Rebound function - calculates the rebound effect for 1, 2, 4 and 6 years
#-------------------------------------------------------------------------

REBOUND <- function(res, nbs, hori) {
  
  #hori1
  re1 <- ((res[[3]][[1]]$irf[1,2] - res[[3]][[1]]$irf[hori[1],2]) / res[[3]][[1]]$irf[1,2])
  rebs <- NULL
  for (ii in 1:nbs) {
    rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[1,2] - res[[3]][[2]][[ii]]$irf[hori[1],2]) / res[[3]][[2]][[ii]]$irf[1,2])
  }
  req1 <- quantile(rebs, prob=c(0.05,0.95))
  
  #hori2
  re2 <- ((res[[3]][[1]]$irf[1,2] - res[[3]][[1]]$irf[hori[2],2]) / res[[3]][[1]]$irf[1,2])
  rebs <- NULL
  for (ii in 1:nbs) {
    rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[1,2] - res[[3]][[2]][[ii]]$irf[hori[2],2]) / res[[3]][[2]][[ii]]$irf[1,2])
  }
  req2 <- quantile(rebs, prob=c(0.05,0.95))
  
  #hori3
  re3 <- ((res[[3]][[1]]$irf[1,2] - res[[3]][[1]]$irf[hori[3],2]) / res[[3]][[1]]$irf[1,2])
  rebs <- NULL
  for (ii in 1:nbs) {
    rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[1,2] - res[[3]][[2]][[ii]]$irf[hori[3],2]) / res[[3]][[2]][[ii]]$irf[1,2])
  }
  req3 <- quantile(rebs, prob=c(0.05,0.95))
  
  #hori4
  re4 <- ((res[[3]][[1]]$irf[1,2] - res[[3]][[1]]$irf[hori[4],2]) / res[[3]][[1]]$irf[1,2])
  rebs <- NULL
  for (ii in 1:nbs) {
    rebs[ii] <- ((res[[3]][[2]][[ii]]$irf[1,2] - res[[3]][[2]][[ii]]$irf[hori[4],2]) / res[[3]][[2]][[ii]]$irf[1,2])
  }
  req4 <- quantile(rebs, prob=c(0.05,0.95))
  
  return(list(hori1 = c(re1, req1), hori2 = c(re2, req2), hori3 = c(re3, req3), hori4 = c(re4, req4)))	
}

#-------
# fastica
#-------
# Create vector with filenames for loop
filename <-"irf.fastica_M.v3.hddcdd.p3.levels.iter1000.RDATA"

# Create table to collect rebound effects
reb_table <- data.frame(matrix(ncol=9, nrow=length(filename)))
colnames(reb_table) <- c("Name", "1 year","0.90 CI", "2 years","0.90 CI", "4 years", "0.90 CI", "6 years", "0.90 CI")

# Create list to collect B matrices
B_list <- list()

# Create list for confidence intervals of B matrix
contemp <- list()

# Loop for the different models
for (i in 1:length(filename)) {
  load(paste(savepath, filename[i], sep=""))
  
  # 1, 2, 4 and 6 years
  if (i < 3) {horiz <- c(4,8,16,24)}else{horiz <- c(12,24,48,72)} 
  
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
  print(plot(res[[3]], lowerq = 0.05, upperq = 0.95))
  dev.off()
  
  # B matrix
  B_list[[i]] <- list(filename[i], round(res[[2]]$B, 2))
  
  # Confidence interval for B matrix
  contemp_sig <- matrix(ncol = ncol(res[[3]]$bootstrap[[1]]$irf), nrow = length(res[[3]]$bootstrap))
  colnames(contemp_sig) <- colnames(res[[3]]$bootstrap[[1]]$irf)
  for (j in 1:length(res[[3]]$bootstrap)) {
    
    contemp_sig[j,] <- unlist(res[[3]]$bootstrap[[j]]$irf[1,])
  }
  
  q10 <- function(x) {quantile(x, probs=c(0.05,0.95))}
  contemp[[i]] <- list(filename[i], apply(contemp_sig, 2, q10))
  
}

# Rebound effects for fastica (Table 3)
reb_table

# Confidence intervals for fastica for B matrix (Table 1 and 2)
contemp


