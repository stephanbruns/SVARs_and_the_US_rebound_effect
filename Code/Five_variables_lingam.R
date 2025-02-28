rm(list=ls())


#------------------------------------
# This code produces the results for 
#	LiNGAM
#	5 variables SVAR (energy-gdp-prices) 
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
require(svars)
require(pbapply)
require(clue)
require(steadyICA)
require(expm)

require(permute)

source(paste(ppath, "Code/sources/sourcedir.R", sep=""))
sourceDir(paste(ppath, "Code/sources/SVARS_mod/", sep=""), FALSE)
sourceDir(paste(ppath, "Code/sources/maincodes/", sep=""), FALSE)
sourceDir(paste(ppath, "Code/sources/lingam_code/code/", sep=""), FALSE)


#-----------------------------------
# Global settings
#-----------------------------------
nbs <- 1000 #bootstrap iterations 
hori_Q <- 80 #horizon quarterly data
hori_M <- 80 #horizon quarterly data
sseed<-42 #for exact reproducibility of the results

loadpath <- paste(ppath, "Data/", sep="") 
savepath <- paste(ppath, "Results/", sep="") 

# function to estimate lingam:
flingam<-function(ures, prune=TRUE, seedl=sseed){
  set.seed(seedl)
  reslg <- lingam(t(ures))
  B0hat <- reslg$Bpruned
  if(prune==FALSE){B0hat <- reslg$Bnopruned}
  k<-ncol(B0hat)
  G0 <- diag(k) - B0hat
  MixM<-solve(G0)
  MixMs<-MixM
  eres<- as.data.frame(t(G0 %*% t(ures)))
  for(i in 1:k){
    MixMs[,i]<-MixM[,i]*sd(eres[,i]) # important: times st.dev.
  }
  ord <- as.vector(reslg$k)
  list(MixMs, ord)
}

# function to check lingam initial conditions stability
f_initcond_stability<-function(ures, ord){
  cord<-matrix(nrow=100, ncol=ncol(ures))
  ck<-cord
  for (i in 1:100){
    set.seed(sseed+i)
    reslg <- lingam(t(ures))
    cord[i,] <- as.vector(reslg$k)
    ck[i,]<-(cord[i,]==ord)
  }
  length(which(apply(ck,1,sum)==ncol(ures)))/100
}

# initialization of the variable iden containing all the identification results
iden<-list(B=NA, A_hat=NA, method="lingam", n=NA, type="const", y=NA, p=NA, K=NA, PIT=FALSE)
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

#-----------------------------------
# 5 variables VAR with e, y, p, i, q
#-----------------------------------
D <- dat
k<-ncol(D)
p<-3

#VAR estimation
VAR.est <- vars::VAR(D, p = p, type="const") # estimation of VAR using OLS
ures<-resid(VAR.est) # reduced form VAR residuals
const<-1:k; for (i in 1:k){const[i]<-coef(VAR.est)[[i]]["const",1]}# constant vector
#roots(VAR.est)
MM<-Acoef(VAR.est)


# Test for normality
print(Gauss_Tests(ures))
normtest::jb.norm.test(ures[,1], nrepl=2000) #u_e residual
normtest::jb.norm.test(ures[,2], nrepl=2000) #u_y residual
normtest::jb.norm.test(ures[,3], nrepl=2000) #u_p residual
normtest::jb.norm.test(ures[,4], nrepl=2000) #u_i residual
normtest::jb.norm.test(ures[,5], nrepl=2000) #u_q residual

### lingam estimation
res<-flingam(ures, seedl=sseed)
Ap<-res[[1]]
Ap[,1]<--Ap[,1]
round(Ap,2) #Mixing matrix of fastICA (confidence intervals are calculated below)
#note:  mean of the bootstrapped mixing matrix calculated below
ord1<-res[[2]] #contemporaneous causal order
ord1 # 4 2 5 1 3 means i-> y -> q -> e -> p 

res<-flingam(ures, seedl=46)
Ap<-res[[1]]
Ap[,1]<--Ap[,1]
round(Ap,2) #Mixing matrix of fastICA (confidence intervals are calculated below)
#note: mean of the bootstrapped mixing matrix calculated below
ord2<-res[[2]] #contemporaneous causal order
ord2 # 5 2 4 1 3 


## check whether lingam is stable across resamling of initial conditions:
f_initcond_stability(ures,ord1) # 0.31 as result means that the causal order 4 2 5 1 3 is 31% stable across resampling of initial conditions 
f_initcond_stability(ures,ord2) # 0.52 means 52% stable of 5 2 4 1 3 

ord<-ord2 #we impose the most stable causal structure
          # for robustness check with the the second most stable structure put ord=ord1
# ord<-ord1 #---robustness check
Pm <-bootstrap_lingam(YY=D, AA=MM, cons=const, ures, niter=1000, sseed=sseed)
Pmo<-Pm[order(Pm[,(ncol(ures)+1)], decreasing=TRUE),]
Pmo # The result 

# 5 2 4 1 3 30% stability
# 2 5 4 1 3 17,5% stability
# 2 4 5 1 3 14,2% stability
# notice that 96.3% of the structure have 2(gdp) ordered before 1(e) and both gdp and e before 3(price) 

#------------------impulse response functions ----------------#
boots<-irf_ci_lingam_mod(YY=D, AA=MM, cons=const, lag_irf=79, ures, niter=1000, ord=ord, sseed=sseed)
str(boots$true)
Apm<-matrix(NA,k,k)
count<-1
for(i in 1:k){
  for (j in 1:k){
    count<-count+1
    Apm[i,j]<-boots$true[[1]][1,count]
  }
}
round(Apm,2) # mixing matrix - bootstrap mean 

#[1,] -1.64 0.48  0.00 0.54 -0.16
#[2,]  0.00 1.22  0.00 0.00  0.22
#[3,]  0.59 1.16 13.66 4.65 -4.65
#[4,]  0.00 1.75  0.00 1.30  0.08
#[5,]  0.00 0.00  0.00 0.00  0.80
#
# confidence intervals
q10 <- function(x) {quantile(x, probs=c(0.05,0.95))}
contemp_sig <- matrix(ncol = ncol(boots$bootstrap[[1]]$irf), nrow = length(boots$bootstrap))
colnames(contemp_sig) <- colnames(boots$bootstrap[[1]]$irf)
for (j in 1:length(boots$bootstrap)) {
contemp_sig[j,] <- unlist(boots$bootstrap[[j]]$irf[1,])
}
c.i.<-apply(contemp_sig, 2, q10)
round(c.i.,2) #### Confidence intervals for mixing matrix

#save data
iden$B<-Apm
iden$A_hat<-matrix(nrow=k, ncol=1+(k*p))
iden$A_hat[,1]<-const
for(i in 1:p){iden$A_hat[,(i*k+2-k):(i*k+1)]<-MM[[i]]}
iden$n<-nrow(D)
iden$y<-as.matrix(t(D))
iden$p<-p; iden$k<-k
res <- list("uniB.lingam_Q.v5.p3.levels", iden, boots)
save(res, file=paste(savepath, "irf.lingam_Q.v5.p3.levels.iter1000.RDATA", sep=""))


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
# 5 variables VAR with e, y, p, i, q 
#----------------------------------
D <- dat[,c("e","y","p","i","q")]
k<-ncol(D)
#VAR estimation
p=3
VAR.est <- vars::VAR(D, p = p, type="const") # estimation of VAR using OLS
ures<-resid(VAR.est) # reduced form VAR residuals
const<-1:k; for (i in 1:k){const[i]<-coef(VAR.est)[[i]]["const",1]}# constant vector
MM<-Acoef(VAR.est)

# Test for normality
print(Gauss_Tests(ures))
normtest::jb.norm.test(ures[,1], nrepl=2000) #u_e residual
normtest::jb.norm.test(ures[,2], nrepl=2000) #u_y residual
normtest::jb.norm.test(ures[,3], nrepl=2000) #u_p residual
normtest::jb.norm.test(ures[,4], nrepl=2000) #u_i residual
normtest::jb.norm.test(ures[,5], nrepl=2000) #u_q residual

### lingam estimation
res<-flingam(ures, seedl=sseed)
Ap<-res[[1]]
Ap[,1]<--Ap[,1]
round(Ap,2) #Mixing matrix of fastICA (confidence intervals are calculated below)
#note:  mean of the bootstrapped mixing matrix calculated below
ord1<-res[[2]] #contemporaneous causal order
ord1 # 5 4 2 1 3 means q-> i -> y -> e -> p 

res<-flingam(ures, seedl=43)
Ap<-res[[1]]
Ap[,1]<--Ap[,1]
round(Ap,2) #Mixing matrix of fastICA (confidence intervals are calculated below)
#note: mean of the bootstrapped mixing matrix calculated below
ord2<-res[[2]] #contemporaneous causal order
ord2 # 5 2 4 1 3 


## check whether lingam is stable across resamling of initial conditions:
f_initcond_stability(ures,ord1) # 0.36 as result means that the causal order 5 4 2 1 3 is 36% stable across resampling of initial conditions 
f_initcond_stability(ures,ord2) # 0.64 means 64% stable of 5 2 4 1 3 (q y i e p)

ord<-ord2 #we impose the most stable causal structure
# for robustness check with the the second most stable structure put ord=ord1
# ord<-ord1 #---robustness check


Pm <-bootstrap_lingam(YY=D, AA=MM, cons=const, ures, niter=1000, sseed=sseed)
Pmo<-Pm[order(Pm[,(ncol(ures)+1)], decreasing=TRUE),]
Pmo # The result 
#[1,] 2 5 4 1 3 353 (y q i e p)
#[2,] 5 4 2 1 3 203 ord 1
#[3,] 4 2 5 1 3 134
#[4,] 5 2 4 1 3 117 ord 2
# notice that 75% of the structure have 2(gdp) ordered before 1(e) and both gdp and e before 3(price) 



#------------------impulse response functions ----------------#
boots<-irf_ci_lingam_mod(YY=D, AA=MM, cons=const, lag_irf=79, ures, niter=1000, ord=ord, sseed=sseed)
str(boots$true)
Apm<-matrix(NA,k,k)
count<-1
for(i in 1:k){
  for (j in 1:k){
    count<-count+1
    Apm[i,j]<-boots$true[[1]][1,count]
  }
}
round(Apm,2) # this is the point estimate of the mix matrix
#
#[,1] [,2] [,3] [,4]  [,5]
#[1,] -1.73 0.01 0.00 0.42 -0.51
#[2,]  0.00 0.64 0.00 0.00  0.04
#[3,] -1.07 0.78 8.63 1.19 -0.14
#[4,]  0.00 0.46 0.00 0.86 -0.03
#[5,]  0.00 0.00 0.00 0.00  0.82
# confidence intervals
q10 <- function(x) {quantile(x, probs=c(0.05,0.95))}
contemp_sig <- matrix(ncol = ncol(boots$bootstrap[[1]]$irf), nrow = length(boots$bootstrap))
colnames(contemp_sig) <- colnames(boots$bootstrap[[1]]$irf)
for (j in 1:length(boots$bootstrap)) {
  contemp_sig[j,] <- unlist(boots$bootstrap[[j]]$irf[1,])
}
c.i.<-apply(contemp_sig, 2, q10)
round(c.i.,2) #### Confidence intervals of Lingam in Table 2

#save data
iden$B<-Apm
iden$A_hat<-matrix(nrow=k, ncol=1+(k*p))
iden$A_hat[,1]<-const
for(i in 1:p){iden$A_hat[,(i*k+2-k):(i*k+1)]<-MM[[i]]}
iden$n<-nrow(D)
iden$y<-as.matrix(t(D))
iden$p<-p; iden$k<-k
res <- list("uniB.lingam_M.v5.p3.levels", iden, boots)
save(res, file=paste(savepath, "irf.lingam_M.v5.p3.levels.iter1000.RDATA", sep=""))

#######---------------------------------------
####### Estimation of rebound effects for both monthly and quarterly data
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
# lingam
#-------
# Create vector with filenames for loop
filename <- c("irf.lingam_Q.v5.p3.levels.iter1000.RDATA",
              "irf.lingam_M.v5.p3.levels.iter1000.RDATA")

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

# Rebound effects for fastica (Table 4)
reb_table

# Confidence intervals for fastica for B matrix 
contemp


