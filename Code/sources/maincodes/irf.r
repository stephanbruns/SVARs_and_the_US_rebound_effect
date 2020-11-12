irf <- function(AA, B0, lag_irf, u_res){
  # NB the first k lines of IRF are the responses of the first k shocks on the first variable
k<-ncol(u_res)
#### structural residuals
   v_res<-u_res
   for(i in 1:nrow(u_res)) {
   v_res[i,]<- (diag(k)-B0) %*% u_res[i,]
   }
###### MA components
  k<-nrow(B0)
  FI<- MA_components(AA,lag_irf)
  PSI<- FI
  for (i in 1:length(FI))  {
  PSI[[i]]<-FI[[i]] %*% solve( diag(k) - B0 )
  }
#### impulse response functions
   IRF<-matrix(nrow=k*k, ncol=lag_irf+1)
   count<-1
   for (i in 1:k){
   for (j in 1:k){
   for (g in 1:(lag_irf+1)){
   IRF[count,g]<- PSI[[g]][i,j]*sd(v_res[,j])
   }
    count<-count+1
   }
   }
IRF
}