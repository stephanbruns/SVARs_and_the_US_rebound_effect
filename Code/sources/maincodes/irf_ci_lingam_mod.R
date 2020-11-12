
irf_ci_lingam_mod <- function(YY, AA, cons, lag_irf, ures, niter, ord, sseed=42, prune=TRUE){

#Data, AA, B0, lag_irf, u_res, 10
  t<-nrow(ures)
  k<-ncol(ures)
  p<-length(AA) # number of lags (AA must be a list)
  
  # initialization of the impulse response function data in list format
  cn<-colnames(ures)
  nt<-lag_irf + 1#number observations impulse response functions
  nvb<-1+k^2
  irfdf<-as.data.frame(matrix(nrow=nt,ncol=nvb))
  irfdf[,1]<-1:nt
  colnames(irfdf)[1]<-"V1"
  count<-1
  for(i in 1:k){
    for(j in 1:k){
      count<-count+1
      colnames(irfdf)[count]<-paste("epsilon[", cn[j], "] %->%", cn[i])
    }}
  basis<-list(true=NA, bootstrap=NA, SE=NA, nboot=nbs, rademacher=NA, point_estimate=NA, boot_mean=NA,
              signrest=NA, sign_complete=NA, sign_part=NA,cov_bs=NA,method="Wild bootstrap")
  basis$true<-list(irf=irfdf)
  basis$bootstrap<-as.list(1:nbs)
  for(i in 1:nbs){
    basis$bootstrap[[i]]<-list(irf=irfdf)
  }
  class(basis)<-"sboot"
  
  
  


MM<-array(0, c(k*k, lag_irf+1, niter))

countinstab<-0
count_iter<-0
maxniter<-niter*10
for (it in 1:maxniter) {
  count_iter<-count_iter+1
  # block bootstrap 
 ll<-4 # length of block
 lb<-t%%ll  #
 mc<- t%/%ll #(t-lb)/ll
 ind<-1: (t-lb)
 set.seed(sseed+it)
 pp<-sample(mc,mc, replace=TRUE)*ll
  for(w in 1:mc){
  ind[(w*ll - ll +1):(w*ll)]<- (pp[w]-ll+1): pp[w]
  }
  if (lb!=0) {
  s1<- 1:lb
  ind<-c(s1, ind + lb)
  }

# wild bootstrap
#ind<-sample(1:t, t, replace=TRUE)
#ind<-1:t
unew<-ures[ind,]
Ynew<-matrix(0,nrow=(t+p),ncol=k)
Ynew[1:p,]<-as.matrix(YY[1:p,])
for(i in (p+1):(t+p)){
for(j in 1:p){
Y<- AA[[j]]%*%t(YY[i-j,])
Ynew[i,]<-Ynew[i,]+Y 
}
Ynew[i,]<- cons + Ynew[i,] + unew[i-p,]
}
Ynew<-as.data.frame(Ynew)

## method of estimation ##
#
# - cointegrating var
#ci <- ca.jo(Ynew, type = "trace", ecdet = "none", K = p)
#cr<-2 ## cointegration rank
#VAR.est <- vec2var(ci, r=cr) 
# AA_new <-VAR.est$A
# unewh<-VAR.est$resid
#
# - unrestricted var
VARspec <- VAR(Ynew,p, type="const") # estimation of VAR using OLS
unewh<-resid(VARspec)
AA_new<-Acoef(VARspec)
#
# identification via lingam
maxit<-25
for (x in 1:maxit){
reslg <- lingam_mute(t(unewh)) # lingam
ord_new <- as.vector(reslg$k)
if (sum(ord == ord_new) ==k){
  print(paste("causal order after iterations", x))
 break 
}
}
if (sum(ord == ord_new) !=k) {
  print("PROBLEM OF STABILITY")
  countinstab<-countinstab+1
  count_iter<-count_iter-1
  } else{
if(prune==TRUE){B0_new<-reslg$Bpruned}
if(prune==FALSE){B0_new<-reslg$Bnopruned}
psi<-irf(AA_new, B0_new, lag_irf, unewh)
#negative energy shock
es<-seq(1,by=k,length=k)
psi[es,]<- - psi[es,]
MM[,,count_iter]<-psi
psit<-as.data.frame(t(psi))  
basis$bootstrap[[count_iter]]$irf[,2:(k^2+1)]<-psit
}
if (count_iter == niter){break}
}
DOWN<-matrix(nrow=k*k, ncol=lag_irf+1)
UP<-matrix(nrow=k*k, ncol=lag_irf+1)
MED<-matrix(nrow=k*k, ncol=lag_irf+1)
for (m in 1:(lag_irf+1)){
for (n in 1:(k*k)){
mg<-MM[n,m,]
DOWN[n,m]<-  quantile(mg,0.05)#quantile(mg,0.025)#quantile(mg,0.005)#quantile(mg,0.025)# #quantile(mg,0.025)# #mean(mg) - 1.96*sd(mg)
MED[n,m]<- mean(mg)#quantile(mg,0.5)
UP[n,m]<- quantile(mg,0.95)#quantile(mg,0.995)#quantile(mg,0.975)  # # quantile(mg,0.025)quantile(mg,0.975) #quantile(mg,0.995) #mean(mg) + 1.96*sd(mg) 
}
}
#ris <-as.list(1:3)
#ris[[1]]<-DOWN
#ris[[2]]<-MED
#ris[[3]]<-UP
#ris
basis$true$irf[,2:(k^2+1)]<-as.data.frame(t(MED))
basis
}