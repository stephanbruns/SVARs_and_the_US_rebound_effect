
bootstrap_lingam<-function(YY, AA, cons, ures, niter, sseed=42){
k<-ncol(ures) 

allP<-rbind(1:k,allPerms(k))

t<-nrow(ures)
k<-ncol(ures)
p<-length(AA) # number of lags (AA must be a list)

OM<-matrix(0,nrow=niter, ncol=k)
for (it in 1:niter) {
  
  # block bootstrap 
  ll<-4#p # length of block
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
# unewh<-VAR.est$resid
#
# - unrestricted var
VARspec <- VAR(Ynew,p, type="const") # estimation of VAR using OLS
unewh<-resid(VARspec)

# lingam
#
  reslg <- lingam_mute(t(unewh))
  ord_new <- as.vector(reslg$k)
OM[it,]<-ord_new
}
  
ct<-rep(0,nrow(allP))
for(i in 1:nrow(allP)){
for (j in 1:niter){
  if(sum(OM[j,]==allP[i,])==k){ct[i]<-ct[i]+1}
}
}
#max(ct)
#allP[which.max(ct),] # maximum stable orderx
cbind(allP, ct)
}
